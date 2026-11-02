"""Local-solver workflow checks on synthetic FITS; no sky data or network."""

from pathlib import Path
import subprocess

from astropy.io import fits
from astropy.wcs import WCS, Sip
from ccdproc import ImageFileCollection
import numpy as np
import pytest

from photozpy import ImageCollectionGroup, PlateSolving
from photozpy.calibration import plate_solving as module


def image(directory, name='master.fits', **cards):
    directory.mkdir(parents=True, exist_ok=True)
    header = fits.Header({
        'IMTYPE': 'Master Light', 'FIELD': 'field_a', 'FILTER': "SDSS_g'",
        'BUNIT': 'ADU', 'RA': '07 13 59.02', 'DEC': '+01 47 50.00',
        'NCOMBINE': 10, 'EXPTIME': 100., 'REGREF': 'g_001.fits',
    })
    header.update(cards)
    path = directory / name
    fits.HDUList([
        fits.PrimaryHDU(np.arange(64, dtype=np.uint16).reshape(8, 8), header),
        fits.ImageHDU(np.eye(8, dtype=np.uint8), name='MASK'),
        fits.ImageHDU(np.fliplr(np.eye(8, dtype=np.uint8)), name='CMASK'),
        fits.ImageHDU(np.ones((2, 2)), name='EXTRA'),
    ]).writeto(path, checksum=True)
    return path


def wcs_header(*, sip=False):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [4.5, 4.5]
    wcs.wcs.crval = [108.49591667, 1.79722222]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.cd = np.array([[-0.0001, 0.], [0., 0.0001]])
    if sip:
        wcs.wcs.ctype = ['RA---TAN-SIP', 'DEC--TAN-SIP']
        a, b = np.zeros((3, 3)), np.zeros((3, 3))
        a[2, 0], b[0, 2] = 1e-5, -1e-5
        wcs.sip = Sip(a, b, None, None, wcs.wcs.crpix)
    return wcs.to_header(relax=True)


@pytest.fixture
def solver(monkeypatch):
    calls = []
    monkeypatch.setattr(module.shutil, 'which', lambda name: '/local/solve-field')

    def run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        assert Path(cmd[1]).name == 'input.fits'
        assert Path(cmd[1]).parent == kwargs['cwd']
        assert '--nsigma' in cmd
        fits.PrimaryHDU(header=wcs_header()).writeto(cmd[cmd.index('--wcs') + 1])
        Path(cmd[cmd.index('--solved') + 1]).write_bytes(b'\x01')
        return subprocess.CompletedProcess(cmd, 0, 'Field solved', '')

    monkeypatch.setattr(module.subprocess, 'run', run)
    return calls


def test_solution_changes_only_wcs_preserves_pixels_hdus_scaling_and_checksums(tmp_path, solver):
    source = image(tmp_path)
    with fits.open(source, mode='update') as hdul:
        hdul[0].header.update(wcs_header(sip=True))
        hdul[0].header['CROTA2'] = 10.
        hdul[0].header['CTYPE1A'] = 'RA---TAN'
        hdul[0].header['PV1_0'] = 1.
    with fits.open(source) as hdul:
        before = [hdu.data.copy() for hdu in hdul]
    result = PlateSolving.add_wcs_locally_for_file(source, overwrite=True)
    assert result.returncode == 0
    cmd, _ = solver[0]
    assert float(cmd[cmd.index('--ra') + 1]) == pytest.approx(108.49591667)
    assert float(cmd[cmd.index('--dec') + 1]) == pytest.approx(1.79722222)
    assert not Path(cmd[1]).parent.exists()
    with fits.open(source, checksum=True) as hdul:
        assert [hdu.name for hdu in hdul] == ['PRIMARY', 'MASK', 'CMASK', 'EXTRA']
        for hdu, original in zip(hdul, before):
            np.testing.assert_array_equal(hdu.data, original)
            assert hdu.verify_checksum() == 1
            assert hdu.verify_datasum() == 1
        header = hdul[0].header
        assert header['BZERO'] == 32768
        assert header['NCOMBINE'] == 10 and header['EXPTIME'] == 100.
        assert not any(key in header for key in ('A_ORDER', 'A_2_0', 'B_ORDER', 'CROTA2', 'CTYPE1A', 'PV1_0'))
        wcs = WCS(header)
        np.testing.assert_allclose(wcs.all_pix2world([[3.5, 3.5]], 0)[0], [108.49591667, 1.79722222])


def test_separate_output_legacy_clean_false_and_explicit_overwrite(tmp_path, solver):
    source = image(tmp_path)
    original = source.read_bytes()
    with pytest.raises(FileExistsError):
        PlateSolving.add_wcs_locally_for_file(source)
    assert not solver
    PlateSolving.add_wcs_locally_for_file(source, clean=False)
    assert (tmp_path / 'mastersolved.fits').is_file()
    assert source.read_bytes() == original


@pytest.mark.parametrize('failure', ['exit', 'unsolved', 'false_marker', 'no_wcs', 'invalid_wcs', 'singular_wcs', 'timeout'])
def test_failure_preserves_input_and_unrelated_neighbor_files(tmp_path, solver, monkeypatch, failure):
    source = image(tmp_path)
    original = source.read_bytes()
    neighbor = tmp_path / 'other.wcs'
    neighbor.write_text('do not delete')
    original_run = module.subprocess.run

    def run(cmd, **kwargs):
        result = original_run(cmd, **kwargs)
        solved, wcs = Path(cmd[cmd.index('--solved') + 1]), Path(cmd[cmd.index('--wcs') + 1])
        if failure == 'exit':
            return subprocess.CompletedProcess(cmd, 1, '', 'bad solver')
        if failure == 'timeout':
            raise subprocess.TimeoutExpired(cmd, 1)
        if failure == 'unsolved':
            solved.unlink()
        elif failure == 'false_marker':
            solved.write_bytes(b'\x00')
        elif failure == 'no_wcs':
            wcs.unlink()
        elif failure == 'invalid_wcs':
            fits.PrimaryHDU().writeto(wcs, overwrite=True)
        elif failure == 'singular_wcs':
            header = wcs_header()
            header['CDELT1'] = 0.
            fits.PrimaryHDU(header=header).writeto(wcs, overwrite=True)
        return result

    monkeypatch.setattr(module.subprocess, 'run', run)
    with pytest.raises((RuntimeError, ValueError)):
        PlateSolving.add_wcs_locally_for_file(source, overwrite=True, timeout=1)
    assert source.read_bytes() == original
    assert neighbor.read_text() == 'do not delete'
    assert set(tmp_path.iterdir()) == {source, neighbor}
    assert not Path(solver[0][0][1]).parent.exists()


def test_sip_solution_retains_distortion(tmp_path, solver, monkeypatch):
    source = image(tmp_path)
    original_run = module.subprocess.run

    def run(cmd, **kwargs):
        result = original_run(cmd, **kwargs)
        fits.PrimaryHDU(header=wcs_header(sip=True)).writeto(cmd[cmd.index('--wcs') + 1], overwrite=True)
        return result

    monkeypatch.setattr(module.subprocess, 'run', run)
    PlateSolving.add_wcs_locally_for_file(source, overwrite=True)
    assert WCS(fits.getheader(source)).sip is not None


def test_solver_time_metadata_does_not_replace_observation_time(tmp_path, solver, monkeypatch):
    source = image(tmp_path, **{'DATE-OBS': '2024-02-16T04:40:45', 'MJD-OBS': 60356.1, 'TIMESYS': 'UTC'})
    original_run = module.subprocess.run

    def run(cmd, **kwargs):
        result = original_run(cmd, **kwargs)
        with fits.open(cmd[cmd.index('--wcs') + 1], mode='update') as hdul:
            hdul[0].header['MJD-OBS'] = 1.
            hdul[0].header['TIMESYS'] = 'TAI'
        return result

    monkeypatch.setattr(module.subprocess, 'run', run)
    PlateSolving.add_wcs_locally_for_file(source, overwrite=True)
    header = fits.getheader(source)
    assert header['DATE-OBS'] == '2024-02-16T04:40:45'
    assert header['MJD-OBS'] == 60356.1
    assert header['TIMESYS'] == 'UTC'


def test_group_batch_forwarding_progress_and_partial_failures(tmp_path, solver, monkeypatch, caplog):
    a = image(tmp_path / 'a', 'one.fits')
    b = image(tmp_path / 'b', 'two.fits')
    image(tmp_path / 'a', 'bias.fits', IMTYPE='Master Bias')
    b_original = b.read_bytes()
    inputs = ImageCollectionGroup((a.parent, b.parent))
    progress_calls = []

    class Progress:
        def __init__(self, values, **kwargs):
            self.values = values
            assert kwargs['desc'] == 'Plate solving'

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def __iter__(self):
            return iter(self.values)

        def set_postfix_str(self, text):
            progress_calls.append(text)

    original_run = module.subprocess.run

    def run(cmd, **kwargs):
        result = original_run(cmd, **kwargs)
        if len(solver) == 2:
            return subprocess.CompletedProcess(cmd, 1, '', 'no match')
        return result

    monkeypatch.setattr(module, 'tqdm', Progress)
    monkeypatch.setattr(module.subprocess, 'run', run)
    plate = PlateSolving(inputs)
    output = tmp_path / 'output'
    result = plate.add_wcs_locally(
        detect_threshold=6, radius=2, extra_args=['--downsample', '2'],
        save_location=output, timeout=45,
    )
    assert isinstance(result, ImageCollectionGroup)
    assert result is plate.image_collection
    assert plate.solved_paths == (output / a.name,)
    assert b in plate.failed_paths
    assert b.read_bytes() == b_original
    assert not (output / b.name).exists()
    assert progress_calls == ['Solving one.fits', 'Solving two.fits']
    assert 'two.fits' in caplog.text
    for cmd, kwargs in solver:
        assert cmd[cmd.index('--nsigma') + 1] == '6.0'
        assert cmd[cmd.index('--radius') + 1] == '2.0'
        assert cmd[-2:] == ['--downsample', '2']
        assert kwargs['timeout'] == 45.
    plate.failed_paths.clear()
    assert b in plate.failed_paths


@pytest.mark.parametrize('extra_args', ['--dir /tmp', '--new-fits=bad.fits', '--di=/tmp', 'other.fits', '--downsample', '-D /tmp', '--ra 123'])
def test_reserved_or_malformed_extra_options_rejected_before_solver(tmp_path, solver, extra_args):
    source = image(tmp_path)
    with pytest.raises(ValueError):
        PlateSolving.add_wcs_locally_for_file(source, extra_args=extra_args, overwrite=True)
    assert not solver


def test_ra_units_and_missing_pointing():
    assert PlateSolving._pointing(fits.Header(), None, None, 'auto') is None
    assert PlateSolving._pointing(fits.Header(), 7.5, -30, 'auto') == (7.5, -30.)
    assert PlateSolving._pointing(fits.Header(), 7.5, -30, 'hourangle')[0] == pytest.approx(112.5)
    with pytest.raises(ValueError, match='both ra and dec'):
        PlateSolving._pointing(fits.Header(), 7.5, None, 'auto')


def test_missing_executable_does_not_prevent_import_or_modify_fits(tmp_path, monkeypatch):
    source = image(tmp_path)
    original = source.read_bytes()
    monkeypatch.setattr(module.shutil, 'which', lambda name: None)
    plate = PlateSolving(ImageFileCollection(tmp_path, keywords='*'), overwrite=True)
    with pytest.raises(FileNotFoundError, match='Local solver not found'):
        plate.add_wcs_locally()
    assert source.read_bytes() == original
    assert not hasattr(plate, 'add_wcs_astrometrynet')


@pytest.mark.parametrize('available', [
    '/path/bin/solve-field', '/opt/homebrew/bin/solve-field', '/usr/local/bin/solve-field',
])
def test_solver_discovery_prefers_path_then_homebrew(monkeypatch, available):
    calls = []

    def which(name):
        calls.append(name)
        if name == 'solve-field' and available == '/path/bin/solve-field':
            return available
        return available if name == available else None

    monkeypatch.setattr(module.shutil, 'which', which)
    assert PlateSolving._executable('solve-field') == available
    assert calls[0] == 'solve-field'
    if available == '/path/bin/solve-field':
        assert calls == ['solve-field']


def test_explicit_missing_solver_does_not_use_another_installation(monkeypatch):
    calls = []

    def which(name):
        calls.append(name)
        return '/opt/homebrew/bin/solve-field' if name == '/opt/homebrew/bin/solve-field' else None

    monkeypatch.setattr(module.shutil, 'which', which)
    with pytest.raises(FileNotFoundError):
        PlateSolving._executable('/custom/missing/solve-field')
    assert calls == ['/custom/missing/solve-field']
