"""Private process boundary for macOS Tk review sessions.

The notebook never creates a Tk root on this path. A short-lived child owns
the native application, so its Dock icon and Tcl state end with the session.
Decision checkpoints also survive a worker failure or notebook interrupt.
"""

from collections.abc import Callable
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def run_isolated(
    images: list[Path],
    rejected_dir: Path,
    adapter_config: dict | None,
    restore_state: Callable[[dict], None],
) -> None:
    """Wait for a GUI child and restore decisions, including on failure."""
    with tempfile.TemporaryDirectory(prefix="photozpy-review-") as directory:
        work = Path(directory)
        request = work / "request.json"
        state = work / "state.json"
        log = work / "worker.log"
        request.write_text(json.dumps({
            "images": [str(path) for path in images],
            "rejected_dir": str(rejected_dir),
            "adapter": adapter_config,
        }), encoding="utf-8")
        environment = os.environ.copy()
        # Use the same checkout as the caller, including editable installs.
        source = str(Path(__file__).resolve().parents[2])
        environment["PYTHONPATH"] = os.pathsep.join(
            part for part in (source, environment.get("PYTHONPATH")) if part
        )
        with log.open("w", encoding="utf-8") as output:
            process = subprocess.Popen(
                [sys.executable, "-m", __name__, str(request), str(state)],
                stdin=subprocess.DEVNULL, stdout=output,
                stderr=subprocess.STDOUT, env=environment,
            )
            try:
                returncode = process.wait()
            finally:
                try:
                    if process.poll() is None:
                        process.terminate()
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait()
                finally:
                    if state.exists():
                        restore_state(json.loads(state.read_text(encoding="utf-8")))
        if returncode != 0:
            details = log.read_text(encoding="utf-8")[-4000:]
            raise RuntimeError(
                f"Review GUI exited with code {returncode}. "
                "Completed decisions are available on the reviewer.\n"
                f"{details}"
            )
        if not state.exists():
            raise RuntimeError("Review GUI exited without returning its decisions.")


def main() -> None:
    """Run exactly one Tk session; never recursively launch another worker."""
    from .adapter import DS9Adapter
    from .image_reviewer import ImageReviewer

    request, state = map(Path, sys.argv[1:])
    config = json.loads(request.read_text(encoding="utf-8"))
    adapter = DS9Adapter(**(config["adapter"] or {}))
    reviewer = ImageReviewer(
        config["images"], config["rejected_dir"],
        adapter=adapter, isolated=False,
    )
    reviewer._state_path = state
    reviewer.run()


if __name__ == "__main__":
    main()
