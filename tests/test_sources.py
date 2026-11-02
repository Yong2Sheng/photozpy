"""Check source/result association without running image processing."""

from dataclasses import replace

from astropy.coordinates import SkyCoord
import pytest

from photozpy.sources import (
    FilterMagnitude,
    ObservationField,
    SourceInfo,
    Sources,
)
from photozpy.telescope import SARA_CT, SARA_RM


@pytest.fixture
def infos():
    coordinate = SkyCoord(150, 2, unit="deg")
    field = ObservationField(
        field_name="science_field",
        telescope=SARA_CT,
        nominal_pointing=coordinate,
    )
    standard_field = replace(field, field_name="standard_field")
    science = SourceInfo(
        source_name="target",
        source_role="science",
        source_coordinate=coordinate,
        observation_field=field,
    )
    standard = SourceInfo(
        source_name="standard",
        source_role="standard",
        source_coordinate=coordinate,
        observation_field=standard_field,
    )
    return science, standard


def measurement(system="instrumental", telescope=SARA_CT):
    return FilterMagnitude(
        telescope=telescope,
        filter_names=("SDSS g",),
        magnitudes=(18.2,),
        errors=(0.03,),
        significance=(12.0,),
        magnitude_system=system,
    )


def test_metadata_preserved_and_all_results_initially_undefined(infos):
    sources = Sources.from_list(iter(infos))
    assert sources.source_infos == infos
    assert sources.source_names == ("target", "standard")
    assert sources.nsources == len(sources) == 2
    for single, info in zip(sources, infos):
        assert isinstance(single, Sources)
        assert len(single) == 1
        assert single.source_info is info
        for kind in ("instrumental", "reference", "calibrated"):
            assert single.get_magnitudes(kind=kind) is None


def test_index_and_slice_keep_python_order_and_bounds(infos):
    sources = Sources(*infos)
    assert sources[-1].source_info is infos[-1]
    assert sources[:].source_infos == infos
    assert sources[::-1].source_infos == infos[::-1]
    assert sources[1:].source_infos == infos[1:]
    assert len(sources[:0]) == 0
    with pytest.raises(IndexError):
        sources[2]
    with pytest.raises(TypeError):
        sources[1.5]


def test_results_stay_with_source_through_selection_and_clearing(infos):
    sources = Sources(*infos)
    selected = sources.by_role("science").for_field("science_field")
    result = measurement()
    selected.set_magnitudes(kind="instrumental", magnitudes=result)
    assert sources["target"].get_magnitudes(kind="instrumental") is result
    assert sources[::-1][-1].get_magnitudes(kind="instrumental") is result
    assert sources["standard"].get_magnitudes(kind="instrumental") is None
    sources.set_magnitudes(
        source_name="target", kind="instrumental", magnitudes=None
    )
    assert selected.get_magnitudes(kind="instrumental") is None


def test_three_results_coexist_without_conversion_or_overwrite(infos):
    source = Sources(infos[0])
    results = {
        "instrumental": measurement(),
        "reference": measurement("Vega"),
        "calibrated": measurement("AB"),
    }
    for kind, result in results.items():
        source.set_magnitudes(kind=kind, magnitudes=result)
    for kind, result in results.items():
        assert source.get_magnitudes(kind=kind) is result
        filter_name = result.telescope.resolve_filter_alias("SDSS g")
        index = result.filter_names.index(filter_name)
        assert result.magnitudes[index].value == 18.2


def test_same_name_in_different_fields_requires_explicit_selection(infos):
    repeated = replace(infos[1], source_name="target")
    sources = Sources(infos[0], repeated)
    with pytest.raises(ValueError, match="ambiguous"):
        sources["target"]
    with pytest.raises(ValueError, match="ambiguous"):
        sources.set_magnitudes(
            source_name="target", kind="instrumental", magnitudes=measurement()
        )
    assert sources.for_field("science_field")["target"].source_info is infos[0]
    result = measurement("AB")
    sources.set_magnitudes(
        source_name="target", field_name="standard_field",
        kind="reference", magnitudes=result,
    )
    assert sources.get_magnitudes(
        source_name="target", field_name="standard_field", kind="reference"
    ) is result
    assert sources[0].get_magnitudes(kind="reference") is None


def test_duplicate_identity_and_legacy_objects_rejected(infos):
    with pytest.raises(ValueError, match="Duplicate"):
        Sources(infos[0], replace(infos[0]))
    with pytest.raises(TypeError, match="SourceInfo"):
        Sources(object())
    with pytest.raises(TypeError, match="from_list"):
        Sources(infos)


def test_empty_and_missing_selections_are_explicit(infos):
    sources = Sources(*infos)
    assert not Sources()
    assert not sources.for_field("missing")
    assert sources.by_role("standard").source_infos == (infos[1],)
    with pytest.raises(KeyError):
        sources["missing"]
    with pytest.raises(KeyError):
        Sources().source_info
    with pytest.raises(ValueError, match="ambiguous"):
        sources.source_info
    with pytest.raises(ValueError, match="source_role"):
        sources.by_role("target")


@pytest.mark.parametrize(
    "kind,system",
    [("instrumental", "AB"), ("reference", "instrumental"),
     ("calibrated", "instrumental")],
)
def test_incompatible_system_rejected_without_replacing_result(
    infos, kind, system
):
    source = Sources(infos[0])
    with pytest.raises(ValueError, match="incompatible"):
        source.set_magnitudes(kind=kind, magnitudes=measurement(system))
    assert source.get_magnitudes(kind=kind) is None


def test_bad_result_type_kind_and_telescope_rejected(infos):
    source = Sources(infos[0])
    with pytest.raises(TypeError, match="FilterMagnitude"):
        source.set_magnitudes(kind="instrumental", magnitudes=[18.2])
    with pytest.raises(ValueError, match="Telescope"):
        source.set_magnitudes(
            kind="instrumental", magnitudes=measurement(telescope=SARA_RM)
        )
    with pytest.raises(ValueError, match="kind"):
        source.set_magnitudes(kind="unknown", magnitudes=None)
    with pytest.raises(ValueError, match="kind"):
        source.get_magnitudes(kind="unknown")
    assert source.get_magnitudes(kind="instrumental") is None


def test_public_exports():
    from photozpy import Sources as PublicSources

    assert PublicSources is Sources
