"""Small DS9/XPA adapter used by interactive review workflows.

DS9 remains responsible for FITS display, scaling, colormap handling, zooming,
and region editing.  This module only provides a thin Python wrapper around
``xpaset`` and ``xpaget`` so higher-level review code can issue DS9 commands
without constructing subprocess calls directly.
"""

from pathlib import Path
import subprocess

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
        Kept as an adapter-level option for future display policy work.  The
        current implementation applies the scale arguments passed to
        :meth:`load_fits` directly.
    """

    def __init__(
        self,
        software_name: str = "ds9",
        reset_scale_on_load: bool = True,
    ) -> None:
        self.software_name: str = software_name
        self.commands: XPACommands = find_xpa_commands()
        self.reset_scale_on_load: bool = reset_scale_on_load

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
        scale_mode: str = "zscale",
        scale_function: str = "asinh",
    ) -> None:
        """Load a FITS image into DS9 and apply default display settings.

        Paths are resolved and quoted before being sent to DS9 because DS9's
        XPA parser treats spaces in unquoted paths as separate command tokens.
        """

        self.set("fits", _quote_path(path))
        self.set_scale_mode(scale_mode)
        self.set_scale_function(scale_function)
        self.zoom_to_fit()

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
