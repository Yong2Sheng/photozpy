# tests/test_version.py
import photozpy


def test_version_import():
    assert isinstance(photozpy.__version__, str)
