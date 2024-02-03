import os
from datetime import datetime
from pathlib import Path
from urllib.parse import urljoin

import pytz
import requests
import spiceypy
from tqdm import tqdm

"""Absolute path to this package's installation directory."""
project_directory = Path(__file__).resolve().parent

"""NAIF remote server kernel directories."""
naif_path = 'https://naif.jpl.nasa.gov/pub/naif/'


class _Kernel:
    """
    Class to hold information about a kernel and provide the option to download
    a local copy.
    """

    def __init__(self, name: str, kernel_type: str, description: str,
                 remote_path: str, base_path: str = naif_path):
        """
        Parameters
        ----------
        name : str
            The filename of the kernel.
        kernel_type : str
            The kernel type, for example, 'spk' or 'pck'.
        description : str
            A description of the kernel for organizational purposes.
        remote_path : str
            The kernel's remote path on the NAIF repository.
        """
        self._name = name
        self._kernel_type = kernel_type
        self._description = description
        self._relative_path = remote_path
        self._remote_path = remote_path
        self._url = urljoin(base_path, remote_path)
        self._local_path = None

    def __str__(self):
        print_strs = ['NASA NAIF SPICE Kernel', f'   Name: {self._name}',
                      f'   Type: {self._kernel_type.upper()}',
                      f'   Description: {self._description}',
                      f'   Remote path: {self._url}',
                      f'   Local path: {self._local_path}']
        return '\n'.join(print_strs)

    def _check_local_directory(self) -> None:
        """
        Check to see if the expected local directory exists; if not, create it.

        Returns
        -------
        None.
        """
        if not self._local_path.exists():
            self._local_path.mkdir(parents=True)

    @staticmethod
    def _get_remote_timestamp(resp) -> float:
        """
        Get remote kernel timestamp.

        Parameters
        ----------
        resp : requests.models.Response
            NAIF server's response to HTTP request.

        Returns
        -------
        timestamp : float
            The UNIX timestamp of the remote file.
        """
        fmt = '%a, %d %b %Y %H:%M:%S %Z'
        timestamp = datetime.strptime(resp.headers['last-modified'], fmt)
        timestamp = timestamp.timestamp()
        return timestamp

    @staticmethod
    def _get_local_timestamp(filename: Path) -> float:
        """
        Get local kernel timestamp.

        Parameters
        ----------
        filename : Path
            Location of local kernel.

        Returns
        -------
        timestamp : float
            The UNIX timestamp of the remote file.
        """
        timestamp = datetime.utcfromtimestamp(os.path.getmtime(filename))
        timestamp = pytz.utc.localize(timestamp).timestamp()
        return timestamp

    def download(self, verbose: bool = False) -> bool:
        """
        Download kernel and preserve original timestamp.

        Parameters
        ----------
        verbose : bool
            Whether or not to print out process information. Looks cool in the
            terminal window, but also potentially useful for other purposes.

        Returns
        -------
        downloaded : bool
            Whether or not the kernel was downloaded from the remote server.
        """
        self._check_local_directory()
        filename = Path(self._local_path, self._name)
        url = urljoin(self._url, self._name)
        resp = requests.get(url, stream=True)
        total = int(resp.headers.get('content-length', 0))
        remote_timestamp = self._get_remote_timestamp(resp)

        download = False

        # first, check to see if the file exists
        if not filename.exists():
            download = True

        # second, check to see if timestamps match
        if filename.exists() and not download:
            local_timestamp = self._get_local_timestamp(filename)
            if local_timestamp != remote_timestamp:
                download = True

        # finally, check to see if the file sizes match
        if filename.exists() and not download:
            filesize = filename.stat().st_size
            if not filesize == total:
                download = True

        chunk_size = 1024
        tqdm_params = dict(total=total, unit='B', unit_scale=True,
                           unit_divisor=chunk_size, desc=f'      Progress')
        if download:
            print(f'   Downloading {self._url}{self._name}')
            with open(filename, 'wb') as file, tqdm(**tqdm_params) as bar:
                for data in resp.iter_content(chunk_size=chunk_size):
                    size = file.write(data)
                    bar.update(size)
            os.utime(filename, (remote_timestamp, remote_timestamp))
        else:
            if verbose:
                print(f'   {self._name} up to date')

        downloaded = download
        return downloaded

    def set_local_path(self, path: str or Path) -> None:
        """
        Set kernel local path.

        Parameters
        ----------
        path : str or Path
            The local path.

        Returns
        -------
        None.
        """
        self._local_path = Path(path, self._relative_path)

    @property
    def local_path(self) -> Path:
        """
        Get kernel full local file path.
        """
        return Path(self._local_path, self._name)

    @property
    def remote_path(self) -> str:
        """
        Get kernel full remote url.
        """
        return self._remote_path

    @property
    def name(self) -> str:
        """
        Get kernel name.
        """
        return self._name


class _SPICEKernels:
    """
    Wrapper class to check for and download a set of SPICE kernels.
    """

    def __init__(self, local_kernel_path: str or Path, kernel_category: str):
        """
        Parameters
        ----------
        local_kernel_path : str or Path
            The location where the local SPICE kernels are stored. Assumes the
            same nested file structure as https://naif.jpl.nasa.gov/pub/naif/.

        kernel_category : str
            The kernel category, e.g., 'generic' or 'JWST'. This is effectively
            the key for the _all_kernels dictionary, but the result will also
            end up being printed in verbose mode.
        """
        self._local_kernel_path = local_kernel_path
        self._kernel_category = kernel_category
        self._kernels, self._verbose = self._set_local_paths(
            _all_kernels[kernel_category])

    def _set_local_paths(
            self, kernels: list[_Kernel]) -> tuple[list[_Kernel], bool]:
        """
        Set local file paths for all kernels in a set.

        Parameters
        ----------
        kernels : list[_Kernel]
            The list of _Kernel objects.

        Returns
        -------
        kernels : list[_Kernel]
            The same list of _Kernel objects, just with the local file paths
            set.
        verbose : bool
            False if every kernel exists locally. True if one or more do not.
        """
        verbose = False
        for kernel in kernels:
            kernel.set_local_path(self._local_kernel_path)
            if not kernel.local_path.exists():
                verbose = True
        return kernels, verbose

    def download(self) -> None:
        """
        Download a local copy of a SPICE kernel, as needed.

        Returns
        -------
        None.
        """
        if self._verbose:
            print_str = f'Checking {self._kernel_category} SPICE kernels...'
            print(print_str)
        downloads = []
        for kernel in self._kernels:
            downloaded = kernel.download(verbose=self._verbose)
            downloads.append(downloaded)

    def furnish(self) -> None:
        """
        Furnish all kernels in the set.

        Returns
        -------
        None.
        """
        for kernel in self._kernels:
            spiceypy.furnsh(str(kernel.local_path))


"""The most generic of the generic kernels."""
_generic_kernels = [
    _Kernel('de440.bsp', kernel_type='spk/planets',
            description='Planet ephemeris.',
            remote_path='generic_kernels/spk/planets/'),
    _Kernel('gm_de440.tpc', kernel_type='pck',
            description='Mass parameters for planets, natural satellites, '
                        'and other major objects.',
            remote_path='generic_kernels/pck/'),
    _Kernel('pck00010.tpc', kernel_type='pck',
            description='Orientation and size/shape data for natural bodies.',
            remote_path='generic_kernels/pck/'),
    _Kernel('pck00011.tpc', kernel_type='pck',
            description='Orientation and size/shape data for natural bodies.',
            remote_path='generic_kernels/pck/'),
    _Kernel('naif0012.tls', kernel_type='lsk',
            description='Leap seconds kernel.',
            remote_path='generic_kernels/lsk/')
]

"""Mars satellite kernels."""
_mars_satellites = [
    _Kernel('mar097.bsp', kernel_type='spk',
            description='Martian satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/')
]

"""Jupiter satellite kernels."""
_jupiter_satellites = [
    _Kernel('jup365.bsp', kernel_type='spk',
            description='Jovian regular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('jup344.bsp', kernel_type='spk',
            description='Jovian irregular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
]

"""Saturn satellite kernels."""
_saturn_satellites = [
    _Kernel('sat441.bsp', kernel_type='spk',
            description='Saturnian regular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('sat415.bsp', kernel_type='spk',
            description='Saturnian inner satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('sat452.bsp', kernel_type='spk',
            description='Saturnian irregular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/')
]

"""Uranus satellite kernels."""
_uranus_satellites = [
    _Kernel('ura111.bsp', kernel_type='spk',
            description='Uranian regular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('ura115.bsp', kernel_type='spk',
            description='Uranian inner satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('ura116.bsp', kernel_type='spk',
            description='Uranian irregular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/')
]

"""Neptune satellite kernels."""
_neptune_satellites = [
    _Kernel('nep095.bsp', kernel_type='spk',
            description='Neptunian regular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
    _Kernel('nep102.bsp', kernel_type='spk',
            description='Neptunian irregular satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/'),
]

"""Pluto satellite kernels."""
_pluto_satellites = [
    _Kernel('plu058.bsp', kernel_type='spk',
            description='Plutonian satellites ephemeris.',
            remote_path='generic_kernels/spk/satellites/')]

"""Combined generic kernels."""
_all_generic_kernels = (_generic_kernels + _mars_satellites +
                        _jupiter_satellites + _saturn_satellites +
                        _uranus_satellites + _neptune_satellites +
                        _pluto_satellites)

"""JWST position kernels."""
_jwst_kernels = [
    _Kernel('jwst_pred.bsp', kernel_type='spk',
            description='JWST predicted ephemeris.',
            remote_path='JWST/kernels/spk/'),
    _Kernel('jwst_rec.bsp', kernel_type='spk',
            description='JWST reconstructed ephemeris.',
            remote_path='JWST/kernels/spk/'),
]

"""HST position kernels."""
_hst_kernels = [
    _Kernel('hst.bsp', kernel_type='spk',
            description='HST ephemeris up to 1 month in the future.',
            remote_path='HST/kernels/spk/'),
]

"""Dictionary of all SPICE kernels."""
_all_kernels = {
    'generic': _all_generic_kernels,
    'JWST': _jwst_kernels,
    'HST': _hst_kernels
}
