"""Exercise the process boundary with disposable, non-GUI child programs."""

import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import Mock

import pytest

from photozpy.interactive_ds9 import _review_process as worker


@pytest.mark.parametrize("exitcode", [0, 7])
def test_child_exit_restores_decisions(tmp_path, monkeypatch, exitcode):
    """A real child process exits; the caller keeps its PID and decisions."""
    real_popen = subprocess.Popen
    requested = []
    state = {"index": 1, "history": ["reject"],
             "rejected": [str(tmp_path / "rejected" / "image.fits")]}

    def start(command, **kwargs):
        requested.append(command)
        script = (
            "import json,pathlib,sys; "
            "request=json.loads(pathlib.Path(sys.argv[1]).read_text()); "
            "assert request['images'][0].endswith('image.fits'); "
            f"pathlib.Path(sys.argv[2]).write_text({json.dumps(state)!r}); "
            "print('worker diagnostic'); "
            f"sys.exit({exitcode})"
        )
        return real_popen(
            [sys.executable, "-c", script, *command[-2:]], **kwargs,
        )

    monkeypatch.setattr(worker.subprocess, "Popen", start)
    restore = Mock()
    if exitcode:
        with pytest.raises(RuntimeError, match="worker diagnostic"):
            worker.run_isolated([tmp_path / "image.fits"], tmp_path, None, restore)
    else:
        worker.run_isolated([tmp_path / "image.fits"], tmp_path, None, restore)
    restore.assert_called_once_with(state)
    assert requested[0][:3] == [
        sys.executable, "-m", "photozpy.interactive_ds9._review_process",
    ]
    assert not Path(requested[0][-1]).parent.exists()


def test_interrupt_stops_worker_and_restores_checkpoint(tmp_path, monkeypatch):
    state = {"index": 1, "history": ["good"], "rejected": []}
    process = Mock()
    process.wait.side_effect = [KeyboardInterrupt, 0]
    process.poll.return_value = None

    def start(command, **kwargs):
        Path(command[-1]).write_text(json.dumps(state))
        return process

    monkeypatch.setattr(worker.subprocess, "Popen", start)
    restore = Mock()
    with pytest.raises(KeyboardInterrupt):
        worker.run_isolated([], tmp_path, None, restore)
    process.terminate.assert_called_once()
    restore.assert_called_once_with(state)


def test_unresponsive_worker_is_killed_on_interrupt(tmp_path, monkeypatch):
    process = Mock()
    process.wait.side_effect = [
        KeyboardInterrupt, subprocess.TimeoutExpired("review", 5), 0,
    ]
    process.poll.return_value = None
    monkeypatch.setattr(worker.subprocess, "Popen", Mock(return_value=process))
    with pytest.raises(KeyboardInterrupt):
        worker.run_isolated([], tmp_path, None, Mock())
    process.kill.assert_called_once()


def test_success_without_decisions_is_not_reported_as_success(tmp_path, monkeypatch):
    process = Mock()
    process.wait.return_value = 0
    process.poll.return_value = 0
    monkeypatch.setattr(worker.subprocess, "Popen", Mock(return_value=process))
    with pytest.raises(RuntimeError, match="without returning"):
        worker.run_isolated([], tmp_path, None, Mock())
