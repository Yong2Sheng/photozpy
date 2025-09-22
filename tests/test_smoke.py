import pytest


@pytest.mark.smoke
def test_import():
    import photozpy

    assert hasattr(photozpy, "__version__")  # or: assert photozpy is not None
