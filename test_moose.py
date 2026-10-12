#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""What can be checked about this package, with and without MOOSE in the room.

Without it: the input files each analysis writes, the CSV readers, and that a
machine with no MOOSE is told what to do. With it -- in this package's image, or
on a workstation with MOOSE, gmsh and build123d -- the example cantilever is
analysed three ways and compared with its closed forms.

Run with `pytest` from this directory, beside an installed `partcad` (for the
sandbox helpers the scripts import; see `conftest.py`).
"""

import math

import pytest

support = pytest.importorskip("cae_support", reason="needs a partcad with wrappers/cae_support.py")

import cfd_moose  # noqa: E402
import dynamic_moose  # noqa: E402
import fea_moose  # noqa: E402
import moose_common as moose  # noqa: E402
import thermal_moose  # noqa: E402


def test_a_machine_with_no_moose_is_told_what_to_install(monkeypatch):
    monkeypatch.delenv(moose.EXECUTABLE_ENV, raising=False)
    monkeypatch.setattr(moose.shutil, "which", lambda _name: None)
    monkeypatch.setattr(moose, "DIRECTORIES", ())
    monkeypatch.setattr(moose.sys, "prefix", "/nonexistent")
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    with pytest.raises(support.ToolMissing) as raised:
        moose.find_moose()
    message = str(raised.value)
    assert "dockerImage:" in message
    assert "conda.software.inl.gov" in message
    assert moose.EXECUTABLE_ENV in message


def test_the_tables_are_read_by_column(tmp_path):
    (tmp_path / "case_out.csv").write_text("time,hottest,coldest\n0,0,0\n1,40,20\n")
    (tmp_path / "case_out_nodes_0001.csv").write_text("id,x,y,z,T\n0,0,0,0,20\n1,0.1,0,0,40\n")
    assert moose.postprocessors(str(tmp_path)) == {"time": 1.0, "hottest": 40.0, "coldest": 20.0}
    assert moose.vector(str(tmp_path), "nodes")[1]["T"] == 40.0


def test_a_load_is_a_traction_along_minus_z():
    text = fea_moose._input({"self_weight": False}, {"fix0": [1]}, {"load0": [2]}, [{"load": 100.0, "area": 100.0}], 2)
    assert "value = -1000000" in text
    assert "boundary = 'fix0'" in text
    assert "order = SECOND" in text
    assert "value = -0" in text  # no gravity


def test_the_modal_problem_tags_the_mass_for_the_eigen_side():
    text = dynamic_moose._input({"density": 7850}, {"fix0": [1]}, 4, 2)
    assert "extra_vector_tags = 'eigen'" in text and "coefficient = -7850" in text
    assert "n_eigen_pairs = 4" in text
    assert "EigenDirichletBC" in text


def test_heat_is_a_flux_and_convection_meets_the_ambient():
    text = thermal_moose._input(
        {"ambient_temperature": 25.0},
        2,
        [("temperature0", {"temperature": 20.0})],
        [("heat0", {"heat": 1.0, "area": 100.0})],
        [("convection0", {"convection": 10.0})],
    )
    assert "value = 20" in text
    assert "value = 10000" in text
    assert "coefficient = 10" in text and "T_infinity = 25" in text


def test_a_flow_through_an_axis_aligned_face_is_held_normal_to_it():
    assert cfd_moose._normal_only({"faces": [{"normal": [0, 0, -1]}]}) == "0 0 NA"
    assert cfd_moose._normal_only({"faces": [{"normal": [0.6, 0.8, 0]}]}) is None
    assert cfd_moose._normal_only({"faces": []}) is None


def test_nothing_held_is_said_rather_than_solved():
    for module in (fea_moose, dynamic_moose):
        result = module.process("/nonexistent.glb", {"boundary": []})
        assert result["success"] is False and "holds nothing" in result["exception"]


def _here():
    try:
        moose.find_moose()
        return True
    except support.ToolMissing:
        return False


needs_moose = pytest.mark.skipif(not _here(), reason="needs MOOSE")


def _cantilever(analysis):
    build123d = pytest.importorskip("build123d")
    pytest.importorskip("gmsh")
    root = {"port": "root", "location": [[-50, 0, 0], [0, 0, 1], 0]}
    tip = {"port": "tip", "location": [[50, 0, 0], [0, 0, 1], 0]}
    boundary = {
        "fea": [dict(root, fix=True, load=0.0), dict(tip, fix=False, load=100.0)],
        "dynamic": [dict(root, fix=True)],
        "thermal": [
            dict(root, temperature=20.0, heat=None, convection=None),
            dict(tip, temperature=None, heat=1.0, convection=None),
        ],
    }[analysis]
    request = {
        "analysis": analysis,
        "wrapped": build123d.Box(100, 10, 10).wrapped,
        "boundary": boundary,
        "youngs_modulus": 2.0e11,
        "poissons_ratio": 0.3,
        "density": 7850,
        "self_weight": False,
        "mesh_size": 0.03,
        "mesh_order": 2,
        "conductivity": 50.0,
        "modes": 2,
    }
    assert support.prepare(request) == []
    return request


@needs_moose
def test_the_cantilever_deflects_as_a_beam(tmp_path):
    result = fea_moose.process(str(tmp_path / "c.glb"), _cantilever("fea"))
    assert result["success"], result.get("exception")
    deflection = float(result["warnings"][0].split("peak displacement ")[1].split()[0])
    assert deflection == pytest.approx(0.2016, rel=0.01)


@needs_moose
def test_the_cantilever_rings_as_a_beam(tmp_path):
    result = dynamic_moose.process(str(tmp_path / "c.glb"), _cantilever("dynamic"))
    assert result["success"], result.get("exception")
    first = float(result["warnings"][0].split("Natural frequencies ")[1].split(",")[0])
    expected = 1.875104**2 / (2 * math.pi * 0.1**2) * math.sqrt(2.0e11 * 0.01**4 / 12 / (7850 * 1e-4))
    assert first == pytest.approx(expected, rel=0.01)


@needs_moose
def test_the_cantilever_warms_as_a_rod(tmp_path):
    result = thermal_moose.process(str(tmp_path / "c.glb"), _cantilever("thermal"))
    assert result["success"], result.get("exception")
    hottest = float(result["warnings"][0].split(" C to ")[1].split(" C")[0])
    assert hottest == pytest.approx(40.0, rel=1e-3)
