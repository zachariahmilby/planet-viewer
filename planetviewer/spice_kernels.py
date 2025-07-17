import base64
import os
import platform
import shutil
import warnings
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import requests
import spiceypy as spice
from astropy.coordinates import EarthLocation, Latitude, Longitude
from astropy.time import Time
from erfa import ErfaWarning
from spiceypy.utils.exceptions import NotFoundError
from spiceypy.utils.exceptions import SpiceNOSUCHFILE
from tqdm import tqdm

from planetviewer.files import _project_directory

u.imperial.enable()
warnings.filterwarnings("ignore", category=ErfaWarning)

pck_kernel = 'pck00010.tpc'
"""Set planetary constants kernel; the new pck00011.tpc kernel makes some major 
changes to longitudes on Mars and Neptune."""

kernel_path = Path(_project_directory, 'anc', 'kernels')
"""Variable to hold local SPICE kernel path."""


def _find_code(body: str) -> int:
    """
    Function to loop through defined Earth ground station integer ID codes and
    find one available for assignment to an observatory.

    Parameters
    ----------
    body : str
        The name of the body on which the observatory is located.

    Returns
    -------

    """
    start = int(f'{spice.bodn2c(body)}000')
    found = False
    while not found:
        try:
            spice.bodc2n(start)
        except NotFoundError:
            found = True
        else:
            start += 1
    return start


class _SPICEKernel:
    """
    Class to hold information about a SPICE kernel and provide the option to
    download a local copy.
    """

    def __init__(self,
                 name: str,
                 local_path: str or Path,
                 remote_path: str = None,
                 object_name: str = None,
                 id_no: int = None,
                 register_body: bool = False):
        """
        Parameters
        ----------
        name : str
            The name of the kernel.
        local_path : str
            The local path of the kernel.
        remote_path : str or None
            The remote path of the kernel.
        object_name : str
            The name of the ephemeris object. Must be specified if this is an
            SPK kernel for a small body downloaded from Horizons.
        id_no : int, optional
            Ephemeris object ID number. Must be specified if this is an SPK
            kernel for a small body downloaded from Horizons.
        register_body : bool, optional
            Whether to register the body's name and ID number. Must be
            specified if this is an SPK kernel for a small body downloaded from
            Horizons.
        """
        self._name = name
        self._local_path = Path(kernel_path, local_path)
        self._remote_path = remote_path
        self._object_name = object_name
        self._id_no = id_no
        self._register_body = register_body

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
        if self._remote_path is None:
            return

        # make local directory if needed
        if not self._local_path.parent.exists():
            self._local_path.parent.mkdir(parents=True)

        # query the NAIF server for the specific kernel
        try:
            resp = requests.get(self.remote_path, stream=True)
            _ = resp.headers['last-modified']
        except KeyError:
            msg = (f"Kernel {self._name} not found at "
                   f"{self.remote_path.replace(f'{self._name}', '')}. "
                   f"It may have been replaced by a newer version.")
            print(msg)
            old_path = f'{self._remote_path}/a_old_versions/{self._name}'
            resp = requests.get(old_path, stream=True)

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
        if self._register_body:
            spice.boddef(self._object_name, self._id_no)
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
        return self._local_path

    @property
    def remote_path(self) -> str | None:
        """
        Get kernel full remote url.
        """
        return self._remote_path


class _EarthObservatory:
    """
    Define an observatory on the surface of the Earth, specified by a
    name or a set of latitude, longitude and altitude coordinates.
    """
    def __init__(self,
                 name: str,
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
        self._code = _find_code('Earth')
        spice.boddef(self._name, self._code)

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
            lines = [line.replace('<__code__>', f'{self._code}')
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
        if platform.system() == 'Darwin':
            os.system(f'xattr -d com.apple.quarantine {pinpoint_binary} >/dev/null 2>&1')
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
            local_path=path)
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
            lines = [line.replace('<__code__>', f'{self._code}')
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
            local_path=path)
        return [kernel]


def _earth_obs_kernel_from_name(name: str,
                                overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Automatically generate an observing site's kernels using coordinates from
    Astropy.

    Parameters
    ----------
    name : str
        The name of the observatory.
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
                                    altitude=location.height)
    spk_kernel = observatory.make_spk_kernel(overwrite=overwrite)
    tk_kernel = observatory.make_tk_kernel(overwrite=overwrite)
    return spk_kernel + tk_kernel


def _earth_obs_kernel_from_coords(
        name: str,
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
        name=name, latitude=latitude, longitude=longitude, altitude=altitude)
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
    data = pd.merge(data0, data1, how='outer',
                    on=['key','kernel_url','local_directory'])
    kernels = []  # noqa
    for i in range(len(data)):
        name = str(Path(data['kernel_url'].iloc[i]).name)
        kernel = _SPICEKernel(name=name,
                              local_path=data['local_directory'].iloc[i],
                              remote_path=data['kernel_url'].iloc[i])
        kernels.append(kernel)
    return kernels


def _get_spacecraft_kernels() -> list[_SPICEKernel]:
    """
    Return a list of spacecraft kernels from NAIF.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """

    data = pd.read_csv(Path(_project_directory, 'anc', 'spacecraft.dat'),
                       delimiter=',')
    kernels = []  # noqa
    for i in range(len(data)):
        name = str(Path(data['kernel_url'].iloc[i]).name)
        kernel = _SPICEKernel(name=name,
                              local_path=data['local_directory'].iloc[i],
                              remote_path=data['kernel_url'].iloc[i])
        kernels.append(kernel)

    # add 'Galileo' as a body
    spice.boddef('Galileo', -77)

    return kernels


def _get_observatory_kernels(
        overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Construct observatory SPK and FK kernels. A good place to find observatory
    abbreviations and coordinates is https://airmass.org/observatories.

    Parameters
    ----------
    overwrite : bool
        Whether or not to overwrite an existing observer of this name if it
        exists. Default is `True`.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """
    # Set some pre-defined observatories
    kernels = []
    data = pd.read_csv(Path(_project_directory, 'anc', 'observatories.dat'),
                       delimiter=',')
    for i in range(len(data)):
        kernel = _earth_obs_kernel_from_coords(
            name=data['key'].iloc[i],
            latitude=Latitude(data['lat'].iloc[i]*u.deg),
            longitude=Longitude(data['lon'].iloc[i]*u.deg),
            altitude=data['alt'].iloc[i]*u.km,
            overwrite=overwrite)
        kernels.extend(kernel)
    return kernels


def _get_earth_custom_kernels(overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Copy custom kernels for Earth mean equator of date and true equator of date
    into kernels directory.

    Parameters
    ----------
    overwrite : bool
        Whether or not to overwrite an existing observer of this name if it
        exists. Default is `True`.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """

    files = sorted(Path(_project_directory, 'anc').glob('*.tf'))
    kernels = []
    for file in files:
        out = Path(_project_directory, 'anc', 'kernels', 'custom', 'fk',
                   file.name)
        if not out.parent.exists():
            out.parent.mkdir(parents=True)
        if (out.exists() and overwrite) | (not out.exists()):
            shutil.copyfile(file, out)
        kernel = _SPICEKernel(name=file.name,
                              local_path=file)
        kernels.append(kernel)
    return kernels


def _get_small_body_kernels() -> list[_SPICEKernel]:
    """
    Return a list of small body kernels previously generated by a user from
    JPL Horizons.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """
    kernels = []
    anc_file = Path(_project_directory, 'anc', 'small_bodies.dat')
    if anc_file.exists():
        data = pd.read_csv(anc_file, sep=',')
        for i in range(len(data)):
            kernel = _SPICEKernel(name=Path(data['filename'].iloc[i]).name,
                                  local_path=data['filename'].iloc[i],
                                  object_name=data['name'].iloc[i],
                                  id_no=data['spk_id'].iloc[i],
                                  register_body=True)
            kernels.append(kernel)
    return kernels


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
    _get_earth_custom_kernels(overwrite=update)


def _furnish_all_kernels() -> None:
    """
    Wrapper function to furnish all package-dependent SPICE kernels, if
    necessary.

    Returns
    -------
    None
        None.
    """
    kernels = (_get_naif_kernels() +
               _get_spacecraft_kernels() +
               _get_observatory_kernels() +
               _get_earth_custom_kernels() +
               _get_small_body_kernels())
    for kernel in kernels:
        try:
            kernel.furnish()
        except SpiceNOSUCHFILE:
            kernel.download()
            kernel.furnish()


def _furnish_kernels_if_not_in_pool():
    count = spice.ktotal('ALL')
    if count == 0:
        _furnish_all_kernels()
    else:
        kernels = (_get_naif_kernels() +
                   _get_spacecraft_kernels() +
                   _get_observatory_kernels() +
                   _get_earth_custom_kernels() +
                   _get_small_body_kernels())
        kernel_names = np.array([str(kernel.local_path) for kernel in kernels])
        furnished_kernels = []
        for i in range(count):
            furnished_kernels.append(spice.kdata(i, 'ALL')[0])
        missing_kernels = np.setxor1d(furnished_kernels, kernel_names)
        for kernel in missing_kernels:
            ind = np.where(kernel_names == kernel)[0]
            kernels[ind].furnish()


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


def get_small_body(name: str = None,
                   spk_id: int = None,
                   record_no: int = None,
                   update: bool = False,
                   start_date: Time = Time('1600-01-01'),
                   end_date: Time = Time('2499-12-31')) -> None:

    """
    Many small bodies (asteroids and comets) do not come in the default SPICE
    ephemeris files. However, they are available from the JPL Horizons
    ephemeris service. By calling this function, you will download and/or load
    the necessary data to perform many calculations with small bodies. However,
    some types of calculations (for example, sub-observer coordinates) may not
    be defined for the body in question.

    Parameters
    ----------
    name : str
        The name of the small body.
    spk_id : int
        The SPICE ephemeris ID code of the small body.
    record_no : int
        If there are multiple SPK files for a particular object, they will each
        have a unique record number. The query to download the SPK file for the
        object will fail, and an excpetion message will print the different
        record numbers. In order to make calculations for this object, you will
        have to specify the record number you want.
    update : bool
        Small body ephemeris calculations are frequently updated. If you want
        to update the local data with the most recent available from Horizons,
        set this to True.
    start_date : Time
        The desired start date for the SPK file. Default is January 1, 1600.
        Some ephemeris objects may need a more restricted time range.
    end_date : Time
        The desired end date for the SPK file. Default is December 31, 2499.
        Some ephemeris objects may need a more restricted time range.

    Returns
    -------
    None
        None.
    """

    url = 'https://ssd.jpl.nasa.gov/api/horizons.api'

    start_time = start_date.strftime('%Y-%b-%d')
    stop_time = end_date.strftime('%Y-%b-%d')

    url += "?format=json&EPHEM_TYPE=SPK&OBJ_DATA=NO"

    # prefer record_id, then spk_id, otherwise search by name
    if record_no is not None:
        url += (f"&COMMAND='{record_no}%3B'"
                f"&START_TIME='{start_time}'"
                f"&STOP_TIME='{stop_time}'")
    elif spk_id is not None:
        url += (f"&COMMAND='DES%3D{spk_id}%3B'"
                f"&START_TIME='{start_time}'"
                f"&STOP_TIME='{stop_time}'")
    else:
        url += (f"&COMMAND='DES%3D{name}%3B'"
                f"&START_TIME='{start_time}'"
                f"&STOP_TIME='{stop_time}'")

    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'spk_file_id' not in data.keys():
            output = np.array(data['result'].split('\n'))
            if'    No matches found.' in output:
                msg = ("No matches found. Try checking the name with the "
                       "Horizons web app at "
                       "https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html. If "
                       "your object name has a parenthetical portion, do not "
                       "include it. For example: 'C/2025 N1 (ATLAS)' should "
                       "just be 'C/2025 N1'.")
            else:
                if spk_id is not None:
                    msg = f'Multiple objects found for SPK-ID {spk_id}.'
                else:
                    msg = f'Multiple objects found for object {name}.'
                msg += (' Specify a record number using the `record_no` kwarg '
                        'from the following list instead:')
                msg += '\n\n'
                ind0 = np.where(output==' Matching small-bodies: ')[0][0] + 2
                ind1 = len(output) - 3
                for line in output[ind0:ind1]:
                    msg += line + '\n'
            raise Exception(msg)
        spk_id = int(data['spk_file_id'])
        if name is None:
            url2 = f'https://ssd-api.jpl.nasa.gov/sbdb.api?spk={spk_id}'
            data2 = requests.get(url2).json()
            name = data2['object']['fullname'].upper()
        if 'spk' in data:
            ext = ''
            if record_no is not None:
                ext = f'_{record_no}'
            spk_filename = f'{spk_id}{ext}.bsp'
            spk_filename = Path(kernel_path, 'custom', 'small_bodies',
                                spk_filename)
            if spk_filename.exists() & (not update):
                pass
            else:
                if not spk_filename.parent.exists():
                    spk_filename.parent.mkdir(parents=True)
                with open(spk_filename, 'wb') as file:
                    file.write(base64.b64decode(data['spk']))

            anc_file = Path(_project_directory, 'anc', 'small_bodies.dat')
            if not anc_file.exists():
                with open(anc_file, 'w') as file:
                    file.write('name,spk_id,filename\n')
                filename = Path(spk_filename).relative_to(kernel_path)
                line = f'{name},{spk_id},{filename}\n'
                with open(anc_file, 'r') as file:
                    if not line in file.readlines():
                        write = True
                    else:
                        write = False
                if write:
                    with open(anc_file, 'a') as file:
                        file.write(line)
            print(f"Small body '{name}' now available as an ephemeris object.")
        else:
            raise Exception(f'No SPK data for {name} available from Horizons.')
