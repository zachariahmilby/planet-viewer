import os
from datetime import datetime
from pathlib import Path
from urllib.parse import urljoin

import astropy.units as u
import pytz
import requests
import spiceypy
from astropy.coordinates import EarthLocation, Latitude, Longitude
from tqdm import tqdm

from planetviewer.files import _project_directory
from planetviewer.files import kernel_path

naif_url = 'https://naif.jpl.nasa.gov/pub/naif/'
"""NAIF remote server kernel directories."""

pck_kernel = 'pck00010.tpc'
"""Set planetary constants kernel; the new pck00011.tpc makes some major 
changes to longitudes on Mars and Neptune."""


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
        local_timestamp = datetime.utcfromtimestamp(
            os.path.getmtime(self.local_path))
        local_timestamp = pytz.utc.localize(local_timestamp).timestamp()
        if remote_timestamp == local_timestamp:
            return True
        else:
            return False

    def download(self) -> None:
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
            if self._compare_timestamps(resp):
                return

        # set TQDM parameters
        total = int(resp.headers.get('content-length', 0))
        chunk_size = 1024
        tqdm_params = dict(total=total, unit='B', unit_scale=True,
                           unit_divisor=chunk_size, desc=f'      Progress')

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
        spiceypy.furnsh(str(self.local_path))

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

    def _make_frame_def(self) -> Path:
        """
        Make observatory input definitions file for use by `pinpoint`.

        Returns
        -------
        out_file : Path
            The path to the definitions file.
        """
        template = Path(_project_directory, 'anc', 'obs_template.def')
        with open(template, 'r') as file:
            lines = file.readlines()
            lines = [line.replace('<__site__>', self._name.upper())
                     for line in lines]
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

    def _make_spk_kernel(self, overwrite: bool = False) -> Path:
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
        out_file = Path(_project_directory, 'anc', 'kernels', 'observatories',
                        f'{self._name.lower()}.bsp')
        pck = Path(kernel_path, 'generic_kernels', 'pck', pck_kernel)
        if (out_file.exists()) & (overwrite is True):
            os.remove(out_file)
        elif (out_file.exists()) & (overwrite is False):
            return out_file
        else:
            in_file = self._make_frame_def()
            if not out_file.parent.exists():
                out_file.parent.mkdir(parents=True)
            cmd = f'{str(pinpoint_binary)} -def {str(in_file)} ' \
                  f'-spk {str(out_file)} -pck {str(pck)}'
            os.system(cmd)
            os.remove(in_file)
        return out_file

    def make_kernel(self, overwrite: bool = False) -> [_SPICEKernel]:
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
        kernel = _SPICEKernel(name=f'{self._name.lower()}.bsp',
                              location=str(path.parent), remote=False)
        code = 399000 + self._code
        spiceypy.boddef(self._name, code)
        return [kernel]


def _earth_obs_kernel_from_name(name: str,
                                code: int) -> list[_SPICEKernel]:
    """
    Automatically generate an observing site kernel using coordinates from
    Astropy.

    Parameters
    ----------
    name : str
        The name of the observatory.
    code : int
        IAU observatory code, available from
        https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html.

    Returns
    -------
    list[_SPICEKernel]
        A list containing the SPK kernel for this observatory.
    """
    location = EarthLocation.of_site(name)
    observatory = _EarthObservatory(
        name=name, latitude=location.lat, longitude=location.lon,
        altitude=location.height, code=code)
    kernel = observatory.make_kernel()
    return kernel


def _earth_obs_kernel_from_coords(name: str,
                                  code: int,
                                  latitude: Latitude,
                                  longitude: Longitude,
                                  altitude: u.Quantity,
                                  overwrite: bool = False
                                  ) -> list[_SPICEKernel]:
    """
    Generate an observing site kernel using user-specified coordinates.

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
        A list containing the SPK kernel for this observatory.
    """
    observatory = _EarthObservatory(
        name=name, latitude=latitude, longitude=longitude, altitude=altitude,
        code=code)
    kernel = observatory.make_kernel(overwrite=overwrite)
    return kernel


def _get_outer_planet_kernels(
        satellite_kernels: list[str]) -> list[_SPICEKernel]:
    """
    Generate a set of SPICE kernels for Mars, Jupiter, Saturn, Uranus, Neptune
    or Pluto with appropriate local paths.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICEKernel objects for Mars, Jupiter, Saturn, Uranus,
        Neptune or Pluto.
    """
    kernels = []
    for satellite_kernel in satellite_kernels:
        kernel = _SPICEKernel(
            name=satellite_kernel, location='generic_kernels/spk/satellites/')
        kernels.append(kernel)
    return kernels


def _get_spacecraft_kernels(kernel_names: list[str],
                            location: str) -> list[_SPICEKernel]:
    """
    Get a list of SPICE kernels for a spacecraft.

    Parameters
    ----------
    kernel_names : list[str]
        A list of kernel names. Probably predicts and reconstructed orbits.
    location : str
        The kernel's location on the NAIF server relative to the directory
        `https://naif.jpl.nasa.gov/pub/naif/`.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICEKernel objects for the spacecraft.
    """
    kernels = []
    for kernel in kernel_names:
        kernels.append(_SPICEKernel(name=kernel, location=location))
    return kernels


common_kernels = [
    _SPICEKernel(name=pck_kernel, location='generic_kernels/pck/'),
    _SPICEKernel(name='gm_de440.tpc', location='generic_kernels/pck/'),
    _SPICEKernel(name='naif0012.tls', location='generic_kernels/lsk/'),
    _SPICEKernel(name='de440.bsp', location='generic_kernels/spk/planets/')]


def _get_uncommon_kernels() -> tuple[dict[str, list[_SPICEKernel]],
                                     list[str], list[str]]:
    """
    Get a dictionary of all the other kernels aside from the common kernels.
    """
    planetary_kernels = {
        'Mercury': [],
        'Venus': [],
        'Earth': [],
        'Mars': _get_outer_planet_kernels(['mar097.bsp']),
        'Jupiter': _get_outer_planet_kernels(['jup344.bsp', 'jup365.bsp']),
        'Saturn': _get_outer_planet_kernels(
            ['sat415.bsp', 'sat441.bsp', 'sat452.bsp']),
        'Uranus': _get_outer_planet_kernels(
            ['ura111.bsp', 'ura115.bsp', 'ura116.bsp']),
        'Neptune': _get_outer_planet_kernels(
            ['nep095.bsp', 'nep097.bsp', 'nep102.bsp']),
        'Pluto': _get_outer_planet_kernels(['plu058.bsp'])
    }

    spacecraft_kernels = {
        'JWST': _get_spacecraft_kernels(['jwst_pred.bsp', 'jwst_rec.bsp'],
                                        location='JWST/kernels/spk/'),
        'HST': _get_spacecraft_kernels(['hst.bsp'],
                                       location='HST/kernels/spk/'),
        'Galileo': _get_spacecraft_kernels(['gll_951120_021126_raj2021.bsp'],
                                           location='GLL/kernels/spk/'),
        'Juno': _get_spacecraft_kernels(
            ['juno_pred_orbit.bsp', 'juno_rec_orbit.bsp'],
            location='JUNO/kernels/spk/'),
    }

    observatory_kernels = {
        'Keck': _earth_obs_kernel_from_coords(
            'Keck', 568, latitude=Latitude('19d 49m 34.7s N'),
            longitude=Longitude('155d 28m 28.0s W'), altitude=4125 * u.m),
        'Palomar': _earth_obs_kernel_from_coords(
            'Palomar', 675, latitude=Latitude('33d 21m 22.6s N'),
            longitude=Longitude('116d 51m 53.4s W'), altitude=1713 * u.m),
        'ALMA': _earth_obs_kernel_from_name('ALMA', 999)
        }

    # add 'Galileo' as a body
    spiceypy.boddef('Galileo', -77)

    spacecraft = list(spacecraft_kernels.keys())
    observatories = list(observatory_kernels.keys())

    all_kernels = {**planetary_kernels, **spacecraft_kernels,
                   **observatory_kernels}

    return all_kernels, spacecraft, observatories


def load_spice_kernels(download: bool = True) -> None:
    """
    Wrapper function to download/update all SPICE kernels and furnish them.
    If you don't want to download/update them, you don't have to.

    Parameters
    ----------
    download : bool
        Whether or not to download/update kernels from NAIF.

    Returns
    -------
    None
        None.
    """
    for kernel in common_kernels:
        if download:
            kernel.download()
        kernel.furnish()
    for kernels in _get_uncommon_kernels()[0].values():
        for kernel in kernels:
            if download:
                try:
                    kernel.download()
                except KeyError:
                    continue
            kernel.furnish()


def get_available_spacecraft() -> list[str]:
    """
    Get a list of all available spacecraft.

    Returns
    -------
    list[str]
        A list of all available spacecraft observatories.
    """
    return _get_uncommon_kernels()[1]


def get_available_observatories() -> list[str]:
    """
    Get a list of all available observatories.

    Returns
    -------
    list[str]
        A list of all available Earth ground-based observatories.
    """
    return _get_uncommon_kernels()[2]
