import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import astropy.units as u
import requests
import spiceypy as spice
from spiceypy.utils.exceptions import SpiceNOSUCHFILE
import pandas as pd
from astropy.coordinates import EarthLocation, Angle, Latitude, Longitude
from tqdm import tqdm

from planetviewer.files import _project_directory

u.imperial.enable()

naif_url = 'https://naif.jpl.nasa.gov/pub/naif/'
"""NAIF remote server kernel directories."""

pck_kernel = 'pck00010.tpc'
"""Set planetary constants kernel; the new pck00011.tpc kernel makes some major 
changes to longitudes on Mars and Neptune."""

kernel_path = Path(_project_directory, 'anc', 'kernels')
"""Variable to hold local SPICE kernel path."""


class _SPICEKernel:
    """
    Class to hold information about a SPICE kernel and provide the option to
    download a local copy.
    """

    def __init__(self,
                 name: str,
                 location: str,
                 remote: bool = True):
        """
        Parameters
        ----------
        name : str
            The name of the kernel.
        location : str
            The kernel's location on the NAIF server relative to the directory
            `https://naif.jpl.nasa.gov/pub/naif/`.
        remote : bool
            Whether or not the kernel is a NAIF kernel which exists on the
            remote server or is one you made yourself.
        """
        self._name = name
        self._local_directory = Path(kernel_path, location)
        self._remote_directory = self._parse_url(location, remote)

    @staticmethod
    def _parse_url(location: str,
                   remote: bool) -> str | None:
        """
        Create a url to the remote file location if it exists.

        Parameters
        ----------
        location : str
            The kernel's location on the NAIF server relative to the directory
            `https://naif.jpl.nasa.gov/pub/naif/`.
        remote : bool
            Whether or not the kernel is a NAIF kernel which exists on the
            remote server or is one you made yourself.

        Returns
        -------
        str | None
            The url or None, as appropriate.
        """
        if remote is False:
            return None
        else:
            return urljoin(naif_url, location)

    def _compare_timestamps(self,
                            resp: requests.models.Response) -> bool:
        """
        Query the NAIF server and compare the timestamp of the local kernel
        copy to the server kernel copy.

        Parameters
        ----------
        resp : requests.models.Response
            Server response from NAIF for the specific kernel.

        Returns
        -------
        bool
            True if the timestamps match, False if they don't or if there is no
            local kernel copy.
        """
        fmt = '%a, %d %b %Y %H:%M:%S %Z'
        remote_timestamp = datetime.strptime(
            resp.headers['last-modified'], fmt).timestamp()
        if not self.local_path.exists():
            return False
        local_timestamp = datetime.fromtimestamp(
            os.path.getmtime(self.local_path), timezone.utc)
        local_timestamp = local_timestamp.timestamp()
        if remote_timestamp == local_timestamp:
            return True
        else:
            return False

    def download(self,
                 overwrite: bool = False) -> None:
        """
        Download a local copy of the kernel. Skips downloading if the local
        copy matches the remote copy.

        Returns
        -------
        None
            None.
        """
        # if this is a custom kernel, skip downloading
        if self._remote_directory is None:
            return

        # make local directory if needed
        if not self._local_directory.exists():
            self._local_directory.mkdir(parents=True)

        # query the NAIF server for the specific kernel
        resp = requests.get(self.remote_path, stream=True)

        # skip if the local copy of the kernel matches the remote copy
        if self.local_path.exists():
            if self._compare_timestamps(resp) and not overwrite:
                return

        # set TQDM parameters
        total = int(resp.headers.get('content-length', 0))
        chunk_size = 1024
        tqdm_params = dict(total=total,
                           unit='B',
                           unit_scale=True,
                           unit_divisor=chunk_size,
                           desc=f'      Progress')

        # download the kernel
        print(f'   Downloading {self.remote_path}')
        with open(self.local_path, 'wb') as file, tqdm(**tqdm_params) as bar:
            for data in resp.iter_content(chunk_size=chunk_size):
                size = file.write(data)
                bar.update(size)

        # set the local timestamp to match the remote timestamp
        fmt = '%a, %d %b %Y %H:%M:%S %Z'
        remote_timestamp = datetime.strptime(
            resp.headers['last-modified'], fmt).timestamp()
        os.utime(self.local_path, (remote_timestamp, remote_timestamp))

    def furnish(self) -> None:
        """
        Furnish the local copy of the kernel.

        Returns
        -------
        None
            None.
        """
        spice.furnsh(str(self.local_path))

    @property
    def name(self) -> str:
        """
        Get kernel name.
        """
        return self._name

    @property
    def local_path(self) -> Path:
        """
        Get kernel full local file path.
        """
        return Path(self._local_directory, self._name)

    @property
    def remote_path(self) -> str | None:
        """
        Get kernel full remote url.
        """
        if self._remote_directory is not None:
            return urljoin(self._remote_directory, self._name)
        else:
            return None


class _EarthObservatory:
    """
    Define an observatory on the surface of the Earth, specified by a
    name or a set of latitude, longitude and altitude coordinates.
    """
    def __init__(self,
                 name: str,
                 code: int,
                 latitude: Latitude,
                 longitude: Longitude,
                 altitude: u.Quantity):
        """
        Parameters
        ----------
        name : str
            The name of the observing site. If the site is one of those known
            to Astropy, then specifying the geodetic coordinates is not
            required.
        code : int
            IAU observatory code, available from
            https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html.
        latitude : Latitude, optional
            Latitude of a user-specified observing site.
        longitude : Longitude, optional
            Longitude of a user-specified observing site.
        altitude : u.Quantity, optional
            Altitude of a user-specified observing site.
        """
        self._name = name
        self._latitude = latitude
        self._longitude = longitude
        self._altitude = altitude.to(u.km)
        self._code = code

    def _make_spk_def(self) -> Path:
        """
        Make observatory input definitions file for use by `pinpoint`.

        Returns
        -------
        out_file : Path
            The path to the definitions file.
        """
        template = Path(_project_directory, 'anc', 'spk_template.def')
        with open(template, 'r') as file:
            lines = file.readlines()
            name = self._name.upper().replace(' ', '_')
            lines = [line.replace('<__site__>', name) for line in lines]
            lines = [line.replace('<__lat__>', f'{self._latitude.deg}')
                     for line in lines]
            lines = [line.replace('<__lon__>',  f'{self._longitude.deg}')
                     for line in lines]
            lines = [line.replace('<__alt__>',  f'{self._altitude.value}')
                     for line in lines]
            lines = [line.replace('<__code__>', f'399{self._code}')
                     for line in lines]
        out_file = Path(_project_directory, 'anc', 'obs.def')
        with open(out_file, 'w') as file:
            file.writelines(lines)
        return out_file

    def _make_spk_kernel(self,
                         overwrite: bool = False) -> Path:
        """
        Make an SPK kernel for the observatory using the `pinpoint` executable.

        Parameters
        ----------
        overwrite : bool, optional
            Whether or not to overwrite an existing kernel. Really only
            necessary if any parameters have changed (like coordinates or valid
            time ranges).

        Returns
        -------
        out_file : Path
            The path to the observatory SPK kernel binary.
        """
        pinpoint_binary = Path(_project_directory, 'anc', 'pinpoint')
        os.chmod(pinpoint_binary, 0o777)  # ensure binary is executable
        out_file = Path(kernel_path, 'observatories', 'spk',
                        f'{self._name.lower().replace(" ", "")}.bsp')
        pck = Path(kernel_path, 'generic_kernels', 'pck', pck_kernel)
        if (out_file.exists()) & (overwrite is True):
            os.remove(out_file)
        elif (out_file.exists()) & (overwrite is False):
            return out_file
        in_file = self._make_spk_def()
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True)
        cmd = f'{str(pinpoint_binary)} -def {str(in_file)} ' \
              f'-spk {str(out_file)} -pck {str(pck)}'
        os.system(cmd)
        os.remove(in_file)
        return out_file

    # noinspection DuplicatedCode
    def make_spk_kernel(self,
                        overwrite: bool = False) -> [_SPICEKernel]:
        """
        Generate a SPICEKernel object for this observatory's SPK kernel.

        Parameters
        ----------
        overwrite : bool, optional
            Whether or not to overwrite an existing kernel. Really only
            necessary if any parameters have changed (like coordinates or valid
            time ranges).

        Returns
        -------
        list[_SPICEKernel]
            A list containing the SPK kernel for this observatory.
        """
        path = self._make_spk_kernel(overwrite=overwrite)
        kernel = _SPICEKernel(
            name=f'{self._name.lower().replace(" ", "")}.bsp',
            location=str(path.parent), remote=False)
        code = int(f'399{self._code:0>3}')
        spice.boddef(self._name, code)
        return [kernel]

    def _make_tk_kernel(self,
                         overwrite: bool = False) -> Path:
        """
        Make an TK kernel for the observatory.

        Parameters
        ----------
        overwrite : bool, optional
            Whether or not to overwrite an existing kernel. Really only
            necessary if any parameters have changed (like coordinates or valid
            time ranges).

        Returns
        -------
        out_file : Path
            The path to the observatory SPK kernel binary.
        """
        template = Path(_project_directory, 'anc', 'tk_template.def')
        with open(template, 'r') as file:
            lines = file.readlines()
            name = self._name.upper().replace(' ', '_')
            lines = [line.replace('<__site__>', name) for line in lines]
            lines = [
                line.replace('<__colat__>', f'{-(90 - self._latitude.deg)}')
                for line in lines]
            lines = [line.replace('<__lon__>', f'{-self._longitude.deg}')
                     for line in lines]
            lines = [line.replace('<__alt__>', f'{self._altitude.value}')
                     for line in lines]
            lines = [line.replace('<__code__>', f'399{self._code}')
                     for line in lines]
        out_file = Path(kernel_path, 'observatories', 'tk',
                        f'{self._name.lower().replace(" ", "")}.tf')
        if (out_file.exists()) & (overwrite is True):
            os.remove(out_file)
        elif (out_file.exists()) & (overwrite is False):
            return out_file
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True)
        with open(out_file, 'w') as file:
            file.writelines(lines)
        return out_file

    # noinspection DuplicatedCode
    def make_tk_kernel(self,
                        overwrite: bool = False) -> [_SPICEKernel]:
        """
        Generate a SPICEKernel object for this observatory's SPK kernel.

        Parameters
        ----------
        overwrite : bool, optional
            Whether or not to overwrite an existing kernel. Really only
            necessary if any parameters have changed (like coordinates or valid
            time ranges).

        Returns
        -------
        list[_SPICEKernel]
            A list containing the SPK kernel for this observatory.
        """
        path = self._make_tk_kernel(overwrite=overwrite)
        kernel = _SPICEKernel(
            name=f'{self._name.lower().replace(" ", "")}.tf',
            location=str(path.parent), remote=False)
        code = int(f'399{self._code:0>3}')
        spice.boddef(self._name, code)
        return [kernel]

def _earth_obs_kernel_from_name(name: str,
                                code: int,
                                overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Automatically generate an observing site's kernels using coordinates from
    Astropy.

    Parameters
    ----------
    name : str
        The name of the observatory.
    code : int
        IAU observatory code, available from
        https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html.
    overwrite : bool, optional
        Whether or not to overwrite an existing kernel. Really only
        necessary if any parameters have changed (like coordinates or valid
        time ranges).

    Returns
    -------
    list[_SPICEKernel]
        A list containing the SPK and FK kernels for this observatory.
    """
    location = EarthLocation.of_site(name)
    observatory = _EarthObservatory(name=name,
                                    latitude=location.lat,
                                    longitude=location.lon,
                                    altitude=location.height,
                                    code=code)
    spk_kernel = observatory.make_spk_kernel(overwrite=overwrite)
    tk_kernel = observatory.make_tk_kernel(overwrite=overwrite)
    return spk_kernel + tk_kernel


def _earth_obs_kernel_from_coords(
        name: str,
        code: int,
        latitude: Latitude,
        longitude: Longitude,
        altitude: u.Quantity,
        overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Generate an observing site's SPK and FK kernels using user-specified
    coordinates.

    Parameters
    ----------
    name : str
        The name of the observatory.
    code : int
        Your choice for an observatory code.
    latitude : Latitude
        Latitude of a user-specified observing site.
    longitude : Longitude
        Longitude of a user-specified observing site.
    altitude : u.Quantity
        Altitude of a user-specified observing site.
    overwrite : bool, optional
        Whether or not to overwrite an existing kernel. Really only
        necessary if any parameters have changed (like coordinates or valid
        time ranges).

    Returns
    -------
    list[_SPICEKernel]
        A list containing the SPK and FK kernels for this observatory.
    """
    observatory = _EarthObservatory(
        name=name, latitude=latitude, longitude=longitude, altitude=altitude,
        code=code)
    spk_kernel = observatory.make_spk_kernel(overwrite=overwrite)
    tk_kernel = observatory.make_tk_kernel(overwrite=overwrite)
    return spk_kernel + tk_kernel


def _get_naif_kernels() -> list[_SPICEKernel]:
    """
    Return a list of planetary and satellite kernels from NAIF. This includes
    the high-precision ITRF93 (Earth) and MEAN_ME (Moon) frames.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """

    path0 = Path(_project_directory, 'anc', 'general_kernels.dat')
    data0 = pd.read_csv(path0, delimiter=',')
    path1 = Path(_project_directory, 'anc', 'planetary_kernels.dat')
    data1 = pd.read_csv(path1, delimiter=',')
    data = pd.merge(data0, data1, how='outer', on=['key','kernel','path'])
    kernels = []
    for i in range(len(data)):
        kernel = _SPICEKernel(name=data['kernel'].iloc[i],
                              location=data['path'].iloc[i])
        kernels.append(kernel)
    return kernels


def _get_spacecraft_kernels() -> list[_SPICEKernel]:
    """
    Return a list of spacecraft kernels from NAIF.
    """

    data = pd.read_csv(Path(_project_directory, 'anc', 'spacecraft.dat'),
                       delimiter=',')
    kernels = []
    for i in range(len(data)):
        kernel = _SPICEKernel(name=data['kernel'].iloc[i],
                              location=data['path'].iloc[i])
        kernels.append(kernel)

    # add 'Galileo' as a body
    spice.boddef('Galileo', -77)

    return kernels


def _get_observatory_kernels(
        overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Construct observatory SPK kernels.

    Parameters
    ----------
    overwrite : bool
        Whether or not to overwrite an existing observer of this name if it
        exists. Default is `True`.
    """

    """Set some observatory sites. Coordinates taken from JPL Horizons."""
    data = pd.read_csv(Path(_project_directory, 'anc', 'observatories.dat'),
                       delimiter=',')
    kernels = []
    for i in range(len(data)):
        kernel = _earth_obs_kernel_from_coords(
            name=data['key'].iloc[i],
            code=data['code'].iloc[i],
            latitude=Latitude(data['lat'].iloc[i]*u.deg),
            longitude=Longitude(data['lon'].iloc[i]*u.deg),
            altitude=data['alt'].iloc[i]*u.km,
            overwrite=overwrite)
        kernels.extend(kernel)

    return kernels


def make_custom_observer(name: str,
                        latitude: Latitude | Angle,
                        longitude: Longitude | Angle,
                        altitude: u.Quantity,
                        code: int = 1000,
                        overwrite: bool = True) -> None:
    """
    Set a custom observer location on Earth's surface.

    Parameters
    ----------
    name : str
        Name of the observer. Should not have any spaces (use underscore "_"
        instead).
    latitude : Latitude
        The geodetic latitude of the observing site.
    longitude : Longitude
        The geodetic longitude of the observing site.
    altitude : u.Quantity
        The altitude of the observing site. Can be in any length units,
        including feet!
    code : int
        Code to use for the observatory. Default is `1000` but if you define
        multiple custom sites these will have to be unique for each site. Best
        to increase from here, so `1001`, `1002`, etc.
    overwrite : bool
        Whether or not to overwrite an existing observer of this name if it
        exists. Default is `True`.

    Returns
    -------
    None
        None.
    """
    kernel = _earth_obs_kernel_from_coords(name,
                                           code,
                                           latitude=Angle(latitude),
                                           longitude=Angle(longitude),
                                           altitude=altitude,
                                           overwrite=overwrite)
    kernel[0].furnish()


def _download_spice_kernels(update: bool = False) -> None:
    """
    Wrapper function to download all SPICE kernels. You must set the local
    directory where you want the kernels downloaded before calling this
    function. The path is stored in the variable `kernel_path`.

    Parameters
    ----------
    update : bool
        Whether or not to overwrite existing kernels. Default is `False`.

    Returns
    -------
    None
        None.
    """
    for kernel in _get_naif_kernels() + _get_spacecraft_kernels():
        kernel.download()
    _get_observatory_kernels(overwrite=update)


def _furnish_spice_kernels() -> None:
    """
    Wrapper function to furnish all package-dependent SPICE kernels.

    Returns
    -------
    None
        None.
    """
    kernels = (_get_naif_kernels() +
               _get_spacecraft_kernels() +
               _get_observatory_kernels())
    for kernel in kernels:
        try:
            kernel.furnish()
        except SpiceNOSUCHFILE:
            kernel.download()
            kernel.furnish()


def update_spice_kernels() -> None:
    """
    Wrapper function to update all SPICE kernels.

    Returns
    -------
    None
        None.
    """
    print('Checking for updates to SPICE kernels...')
    _download_spice_kernels(update=True)
