import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import astropy.units as u
import requests
import spiceypy
from astropy.coordinates import EarthLocation, Angle, Latitude, Longitude
from tqdm import tqdm

from planetviewer.files import _project_directory

u.imperial.enable()

naif_url = 'https://naif.jpl.nasa.gov/pub/naif/'
"""NAIF remote server kernel directories."""

pck_kernel = 'pck00010.tpc'
"""Set planetary constants kernel; the new pck00011.tpc kernel makes some major 
changes to longitudes on Mars and Neptune."""

kernel_path = ''
"""Variable to hold local SPICE kernel path."""


def set_kernel_path(path: str or Path) -> None:
    """
    Set local kernel path.

    Parameters
    ----------
    path

    Returns
    -------

    """
    global kernel_path
    if kernel_path == '':
        kernel_path = path


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
        out_file = Path(kernel_path, 'observatories',
                        f'{self._name.lower()}.bsp')
        pck = Path(kernel_path, 'generic_kernels', 'pck', pck_kernel)
        if (out_file.exists()) & (overwrite is True):
            os.remove(out_file)
        elif (out_file.exists()) & (overwrite is False):
            return out_file
        in_file = self._make_frame_def()
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True)
        cmd = f'{str(pinpoint_binary)} -def {str(in_file)} ' \
              f'-spk {str(out_file)} -pck {str(pck)}'
        os.system(cmd)
        os.remove(in_file)
        return out_file

    def make_kernel(self,
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
        kernel = _SPICEKernel(name=f'{self._name.lower()}.bsp',
                              location=str(path.parent), remote=False)
        code = int(f'399{self._code:0>3}')
        spiceypy.boddef(self._name, code)
        return [kernel]


def _earth_obs_kernel_from_name(name: str,
                                code: int,
                                overwrite: bool = False) -> list[_SPICEKernel]:
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
    overwrite : bool, optional
        Whether or not to overwrite an existing kernel. Really only
        necessary if any parameters have changed (like coordinates or valid
        time ranges).

    Returns
    -------
    list[_SPICEKernel]
        A list containing the SPK kernel for this observatory.
    """
    location = EarthLocation.of_site(name)
    observatory = _EarthObservatory(name=name,
                                    latitude=location.lat,
                                    longitude=location.lon,
                                    altitude=location.height,
                                    code=code)
    kernel = observatory.make_kernel(overwrite=overwrite)
    return kernel


def _earth_obs_kernel_from_coords(
        name: str,
        code: int,
        latitude: Latitude,
        longitude: Longitude,
        altitude: u.Quantity,
        overwrite: bool = False) -> list[_SPICEKernel]:
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


def _get_naif_kernels() -> list[_SPICEKernel]:
    """
    Return a list of planetary and satellite kernels from NAIF. This includes
    the high-precision ITRF93 (Earth) and MEAN_ME (Moon) frames.

    Returns
    -------
    list[_SPICEKernel]
        A list of SPICE kernels.
    """

    planetary_kernels = {
        'Mercury': [],
        'Venus': [],
        'Earth': [_SPICEKernel(name='earth_1962_240827_2124_combined.bpc',
                               location='generic_kernels/pck/'),
                  _SPICEKernel(name='earth_assoc_itrf93.tf',
                               location='generic_kernels/fk/planets/'),
                  _SPICEKernel(name='moon_pa_de440_200625.bpc',
                               location='generic_kernels/pck/'),
                  _SPICEKernel(name='moon_de440_220930.tf',
                               location='generic_kernels/fk/satellites/')],
        'Mars': [_SPICEKernel(name='mar097.bsp',
                              location='generic_kernels/spk/satellites/')],
        'Jupiter': [_SPICEKernel(name='jup344.bsp',
                                 location='generic_kernels/spk/satellites/'),
                    _SPICEKernel(name='jup346.bsp',
                                 location='generic_kernels/spk/satellites/'),
                    _SPICEKernel(name='jup365.bsp',
                                 location='generic_kernels/spk/satellites/')],
        'Saturn': [_SPICEKernel(name='sat415.bsp',
                                location='generic_kernels/spk/satellites/'),
                   _SPICEKernel(name='sat441.bsp',
                                location='generic_kernels/spk/satellites/'),
                   _SPICEKernel(name='sat454.bsp',
                                location='generic_kernels/spk/satellites/')],
        'Uranus': [_SPICEKernel(name='ura115.bsp',
                                location='generic_kernels/spk/satellites/'),
                   _SPICEKernel(name='ura116.bsp',
                                location='generic_kernels/spk/satellites/'),
                   _SPICEKernel(name='ura117.bsp',
                                location='generic_kernels/spk/satellites/'),
                   _SPICEKernel(name='ura182.bsp',
                                location='generic_kernels/spk/satellites/')],
        'Neptune': [_SPICEKernel(name='nep095.bsp',
                                 location='generic_kernels/spk/satellites/'),
                    _SPICEKernel(name='nep104.bsp',
                                 location='generic_kernels/spk/satellites/'),
                    _SPICEKernel(name='nep105.bsp',
                                 location='generic_kernels/spk/satellites/')],
        'Pluto': [_SPICEKernel(name='plu060.bsp',
                               location='generic_kernels/spk/satellites/')]
    }

    kernels = [
        _SPICEKernel(name=pck_kernel, location='generic_kernels/pck/'),
        _SPICEKernel(name='gm_de440.tpc', location='generic_kernels/pck/'),
        _SPICEKernel(name='naif0012.tls', location='generic_kernels/lsk/'),
        _SPICEKernel(name='de440.bsp', location='generic_kernels/spk/planets/')
    ]

    for key, val in planetary_kernels.items():
        kernels.extend(val)

    return kernels


def _get_spacecraft_kernels() -> list[_SPICEKernel]:
    """
    Return a list of spacecraft kernels from NAIF.
    """

    spacecraft_kernels = {
        'JWST': [_SPICEKernel(name='jwst_pred.bsp',
                              location='JWST/kernels/spk/'),
                 _SPICEKernel(name='jwst_rec.bsp',
                              location='JWST/kernels/spk/')],
        'HST': [_SPICEKernel(name='hst.bsp',
                             location='HST/kernels/spk/')],
        'Galileo': [_SPICEKernel(name='gll_951120_021126_raj2021.bsp',
                                 location='GLL/kernels/spk/')],
        'Juno': [_SPICEKernel(name='juno_pred_orbit.bsp',
                              location='JUNO/kernels/spk/'),
                 _SPICEKernel(name='juno_rec_orbit.bsp',
                              location='JUNO/kernels/spk/')],
        'Cassini': [_SPICEKernel(name='171215R_SCPSEops_97288_17258.bsp',
                                 location='CASSINI/kernels/spk/')],
        'Voyager 1': [_SPICEKernel(name='vgr1.x2100.bsp',
                                   location='VOYAGER/kernels/spk/')],
        'Voyager 2': [_SPICEKernel(name='vgr2.x2100.bsp',
                                   location='VOYAGER/kernels/spk/')]
    }

    # add 'Galileo' as a body
    spiceypy.boddef('Galileo', -77)

    kernels = []
    for key, val in spacecraft_kernels.items():
        kernels.extend(val)

    return kernels


def _get_observatory_kernels(
        overwrite: bool = False) -> list[_SPICEKernel]:
    """
    Get a dictionary of all the other NAIF kernels aside from the common
    kernels. This includes the high-precision ITRF93 (Earth) and MEAN_ME (Moon)
    frames, satellites and select spacecraft.

    Parameters
    ----------
    overwrite : bool
        Whether or not to overwrite an existing observer of this name if it
        exists. Default is `True`.
    """

    """Set some observatory sites. Coordinates taken from JPL Horizons."""
    observatory_kernels = {
        'ALMA': _earth_obs_kernel_from_coords(
            'ALMA', 999,
            latitude=Latitude(-23.029211*u.deg),
            longitude=Longitude(292.2452521*u.deg),
            altitude=5.07489*u.km,
            overwrite=overwrite),
        'Apache_Point': _earth_obs_kernel_from_coords(
            'Apache_Point', 705,
            latitude=Latitude(32.7805613*u.deg),
            longitude=Longitude(254.1794*u.deg),
            altitude=2.79488*u.km,
            overwrite=overwrite),
        'Cahill_Roof': _earth_obs_kernel_from_coords(
            'Cahill_Roof', 997,
            latitude=Latitude(34.1354666 * u.deg),
            longitude=Longitude(241.8734532 * u.deg),
            altitude=800*u.imperial.ft,
            overwrite=overwrite),
        'Keck': _earth_obs_kernel_from_coords(
            'Keck', 568,
            latitude=Latitude(19.8260847*u.deg),
            longitude=Longitude(204.5278*u.deg),
            altitude=4.21024*u.km,
            overwrite=overwrite),
        'Kitt_Peak': _earth_obs_kernel_from_coords(
            'Kitt_Peak', 695,
            latitude=Latitude(32.7805613*u.deg),
            longitude=Longitude(254.1794*u.deg),
            altitude=2.79488*u.km,
            overwrite=overwrite),
        'Lowell': _earth_obs_kernel_from_coords(
            'Lowell', 690,
            latitude=Latitude(35.2018865*u.deg),
            longitude=Longitude(248.3367*u.deg),
            altitude=2.22399*u.km,
            overwrite=overwrite),
        'Maunakea': _earth_obs_kernel_from_coords(
            'Maunakea', 568,
            latitude=Latitude(19.8260847 * u.deg),
            longitude=Longitude(204.5278 * u.deg),
            altitude=4.21024 * u.km,
            overwrite=overwrite),
        'McDonald': _earth_obs_kernel_from_coords(
            'McDonald', 711,
            latitude=Latitude(30.6715043*u.deg),
            longitude=Longitude(255.9785*u.deg),
            altitude=2.10673*u.km,
            overwrite=overwrite),
        'Palomar': _earth_obs_kernel_from_coords(
            'Palomar', 675,
            latitude=Latitude(33.354136*u.deg),
            longitude=Longitude(243.1375*u.deg),
            altitude=1.70489*u.km,
            overwrite=overwrite),
        'Sommers_Bausch': _earth_obs_kernel_from_coords(
            'Sommers_Bausch', 463,
            latitude=Latitude(40.0040628*u.deg),
            longitude=Longitude(254.7375*u.deg),
            altitude=1.66163*u.km,
            overwrite=overwrite),
        'VLT': _earth_obs_kernel_from_coords(
            'VLT', 998,
            latitude=Latitude(-24.6254085*u.deg),
            longitude=Longitude(289.5971949*u.deg),
            altitude=2.635*u.km,
            overwrite=overwrite)
        }

    kernels = []
    for key, val in observatory_kernels.items():
        kernels.extend(val)

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
                                           latitude=latitude,
                                           longitude=longitude,
                                           altitude=altitude,
                                           overwrite=overwrite)
    kernel[0].furnish()


def download_spice_kernels(overwrite: bool = False) -> None:
    """
    Wrapper function to download all SPICE kernels. You must set the local
    directory where you want the kernels downloaded before calling this
    function. The path is stored in the variable `kernel_path`.

    Parameters
    ----------
    overwrite : bool
        Whether or not to overwrite existing kernels. Default is `False`.

    Returns
    -------
    None
        None.
    """

    if kernel_path == '':
        raise RuntimeError('You must set the local kernel path.')

    for kernel in _get_naif_kernels() + _get_spacecraft_kernels():
        kernel.download(overwrite=overwrite)
    _get_observatory_kernels(overwrite=overwrite)


def furnish_spice_kernels() -> None:
    """
    Wrapper function to furnish all package-dependent SPICE kernels.

    Returns
    -------
    None
        None.
    """
    if kernel_path == '':
        raise RuntimeError('You must set the local kernel path.')

    kernels = _get_naif_kernels() + _get_spacecraft_kernels() + _get_observatory_kernels()
    for kernel in kernels:
        kernel.furnish()
