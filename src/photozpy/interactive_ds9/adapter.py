"""Small DS9/XPA adapter used by interactive review workflows.

DS9 remains responsible for FITS display, scaling, colormap handling, zooming,
and region editing.  This module only provides a thin Python wrapper around
``xpaset`` and ``xpaget`` so higher-level review code can issue DS9 commands
without constructing subprocess calls directly.
"""

from pathlib import Path
import math
from shutil import which
import subprocess
import sys
import time

from .dependencies import XPACommands, find_xpa_commands


class DS9CommunicationError(RuntimeError):
    """Raised when an XPA command fails while communicating with DS9."""


class DS9Adapter:
    """Control a running DS9 instance through XPA command line tools.

    Parameters
    ----------
    software_name
        XPA access point name.  The default ``"ds9"`` targets the usual DS9
        instance.
    reset_scale_on_load
        If True, reapply the default scale settings on every image load.
        Default is False: initialize the first image with zscale/asinh and
        leave subsequent manual DS9 scale choices unchanged. Explicit scale
        arguments to load_fits always take precedence.
    """

    def __init__(
        self,
        software_name: str = "ds9",
        reset_scale_on_load: bool = False,
    ) -> None:
        self.software_name: str = software_name
        self.commands: XPACommands = find_xpa_commands()
        self.reset_scale_on_load: bool = reset_scale_on_load
        self._has_loaded_image: bool = False
        self._started_process: subprocess.Popen | None = None

    def verify_image_loaded(self) -> None:
        """Require an actual FITS image, not just successful command delivery.

        Controllers use this on initial display to avoid enabling navigation
        over empty frames after a silently ignored or unsupported load request.
        """
        loaded = self.get('frame', 'has', 'fits').strip().lower()
        if loaded not in ('1', 'yes', 'true'):
            raise DS9CommunicationError('DS9 has no FITS image in the active frame; loading failed.')

    def ensure_running(self, timeout: float = 15.) -> None:
        """Connect to this XPA target, or launch DS9 and await readiness.

        Used by quality review and automatic comparison previews.
        Find ds9 on PATH or in the standard macOS application bundle. Leave
        the DS9 application open when the comparison controls are closed.
        """
        def ready() -> bool:
            try:
                version = self.get('version')
                return bool(version) and 'XPA$ERROR' not in version
            except DS9CommunicationError:
                return False

        if ready():
            return

        executable = which('ds9')
        if executable is None and sys.platform == 'darwin':
            app_binary = Path('/Applications/SAOImageDS9.app/Contents/MacOS/ds9')
            if app_binary.is_file():
                executable = str(app_binary)
        if executable is None:
            raise FileNotFoundError('DS9 executable not found; install DS9 or add ds9 to PATH.')

        self._started_process = subprocess.Popen(
            [executable, '-title', self.software_name],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, start_new_session=True,
        )
        self._has_loaded_image = False
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if ready():
                return
            if self._started_process.poll() is not None:
                break
            time.sleep(0.1)

        raise DS9CommunicationError('DS9 did not become available; open DS9 and retry preview.')

    def set(self, *args: str) -> None:
        """Send an ``xpaset -p`` command to DS9.

        Parameters
        ----------
        *args
            DS9 command tokens passed after the XPA target name.

        Raises
        ------
        DS9CommunicationError
            If ``xpaset`` exits with a non-zero return code.
        """

        command = [self.commands.xpaset, "-p", self.software_name, *args]
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
        )
        if result.returncode != 0:
            raise DS9CommunicationError(result.stderr.strip())

    def get(self, *args: str) -> str:
        """Send an ``xpaget`` command to DS9 and return stripped stdout.

        Parameters
        ----------
        *args
            DS9 query tokens passed after the XPA target name.

        Raises
        ------
        DS9CommunicationError
            If ``xpaget`` exits with a non-zero return code.
        """

        command = [self.commands.xpaget, self.software_name, *args]
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
        )
        if result.returncode != 0:
            raise DS9CommunicationError(result.stderr.strip())
        return result.stdout.strip()

    def load_fits(
        self,
        path: str | Path,
        scale_mode: str | None = None,
        scale_function: str | None = None,
        *,
        configure_display: bool = True,
        hdu: str | int | None = None,
        preserve_view: bool = False,
    ) -> None:
        """Load a FITS image without repeatedly overriding manual scale choices.

        Parameters
        ----------
        path : str or Path
            FITS image to display in the active DS9 frame.
        scale_mode : str or None, optional
            Explicit scale-limit mode, such as zscale or minmax. None uses
            zscale on the first load or when reset_scale_on_load is True;
            otherwise no scale-mode command is sent.
        scale_function : str or None, optional
            Explicit stretch function, such as log or asinh. None uses asinh
            on the first load or when reset_scale_on_load is True; otherwise
            no stretch command is sent, preserving the current DS9 choice.
        configure_display : bool, optional
            Default True retains the normal scale initialization and zoom-to-fit.
            False sends only the image-load command, leaving display geometry
            and scaling to DS9/the comparison controller. Explicit scale
            arguments cannot be combined with False.
        hdu : str, int, or None, optional
            Select a FITS extension by name or index, e.g. 'CMASK'. None uses
            the normal DS9 default image selection.
        preserve_view : bool, optional
            Default False retains zoom-to-fit. True captures the active image's
            zoom and image-coordinate pan before loading, then restores both.
            Requires a loaded image and configure_display=True. Intended for
            navigating registered images of the same target or field.

        Paths are resolved and quoted before being sent to DS9 because DS9's
        XPA parser treats spaces in unquoted paths as separate command tokens.
        Zoom-to-fit is applied unless preserve_view=True. A failed load does not
        mark initialization complete, allowing Retry to initialize the display.
        """

        if not configure_display and (scale_mode is not None or scale_function is not None):
            raise ValueError('Explicit scale settings require configure_display=True.')
        if preserve_view and not configure_display:
            raise ValueError('preserve_view requires configure_display=True.')

        view: tuple[list[str], list[str]] | None = None
        if preserve_view:
            zoom = self.get("zoom").split()
            pan = self.get("pan", "image").split()
            try:
                if len(zoom) not in (1, 2) or len(pan) != 2:
                    raise ValueError("Unexpected zoom/pan dimensions.")
                zoom_values = [float(value) for value in zoom]
                pan_values = [float(value) for value in pan]
                if not all(math.isfinite(value) for value in zoom_values + pan_values):
                    raise ValueError("Non-finite zoom/pan values.")
                if any(value <= 0 for value in zoom_values):
                    raise ValueError("Zoom must be positive.")
            except ValueError as error:
                raise DS9CommunicationError(f"Cannot capture DS9 view: {error}") from error
            view = (zoom, pan)

        filename = _quote_path(path)
        if hdu is not None:
            valid_index = isinstance(hdu, int) and not isinstance(hdu, bool) and hdu >= 0
            valid_name = isinstance(hdu, str) and bool(hdu) and all(
                char.isascii() and (char.isalnum() or char in '_-') for char in hdu
            )
            if not (valid_index or valid_name):
                raise ValueError('hdu must be a nonnegative index or a simple FITS extension name.')
            filename = filename[:-1] + f'[{hdu}]"'

        self.set('fits', filename)
        if not configure_display:
            return

        initialize_scale: bool = (
            self.reset_scale_on_load or not self._has_loaded_image
        )
        if scale_mode is not None or initialize_scale:
            self.set_scale_mode("zscale" if scale_mode is None else scale_mode)

        if scale_function is not None or initialize_scale:
            self.set_scale_function(
                "asinh" if scale_function is None else scale_function
            )

        if view is None:
            self.zoom_to_fit()
        else:
            zoom, pan = view
            self.set("zoom", "to", *zoom)
            self.set("pan", "to", *pan, "image")
        self._has_loaded_image = True

    def set_zscale(self) -> None:
        """Set DS9's scale mode to ``zscale``."""

        self.set("scale", "mode", "zscale")

    def set_cmap(self, cmap: str) -> None:
        """Set the active DS9 colormap by name."""

        self.set("cmap", cmap)

    def load_region(self, path: str | Path) -> None:
        """Load a DS9 region file into the active frame."""

        region_path = Path(path).expanduser().resolve()
        self.set("region", "load", _quote_path(region_path))

    def add_regions(self, region_text: str) -> None:
        """Send DS9-format regions through XPA stdin, without a temporary file.

        The caller supplies the coordinate-system declaration and region
        properties. Existing regions remain until clear_regions is called.
        """
        result = subprocess.run(
            [self.commands.xpaset, self.software_name, "regions"],
            input=region_text,
            check=False,
            text=True,
            capture_output=True,
        )
        if result.returncode != 0:
            raise DS9CommunicationError(result.stderr.strip())

    def save_region(self, path: str | Path) -> None:
        """Save the active DS9 regions to a region file."""

        region_path = Path(path).expanduser().resolve()
        self.set("region", "save", _quote_path(region_path))

    def clear_regions(self) -> None:
        """Delete all regions from the active DS9 frame."""

        self.set("region", "delete")

    def set_scale_mode(self, mode: str) -> None:
        """Set DS9's scale limit mode, for example ``zscale`` or ``minmax``."""

        self.set("scale", "mode", mode)

    def set_scale_function(self, function: str) -> None:
        """Set DS9's stretch function, for example ``linear`` or ``asinh``."""

        self.set("scale", function)

    def set_scale_limits(self, lower: float, upper: float) -> None:
        """Set explicit lower and upper DS9 scale limits."""

        self.set("scale", "limits", str(lower), str(upper))

    def set_cmap_parameters(
        self,
        contrast: float = 1.0,
        bias: float = 0.5,
    ) -> None:
        """Set DS9 colormap contrast and bias.

        The values correspond to DS9's colormap controls, not to matplotlib
        normalization objects.  They are useful for display tuning only and do
        not modify FITS data.
        """

        if not 0.0 <= contrast <= 10.0:
            raise ValueError("contrast must be between 0 and 10")
        if not 0.0 <= bias <= 1.0:
            raise ValueError("bias must be between 0 and 1")

        self.set("cmap", str(contrast), str(bias))

    def zoom_to_fit(self) -> None:
        """Ask DS9 to fit the current image into the display window."""

        self.set("zoom", "to", "fit")


def _quote_path(path: str | Path) -> str:
    """Return a DS9-safe quoted absolute path string.

    ``subprocess.run(..., shell=False)`` already protects the local shell side.
    This quoting is for DS9's own XPA command parser, which still needs paths
    containing spaces to arrive as one quoted token.
    """

    resolved = Path(path).expanduser().resolve()
    escaped = str(resolved).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'
