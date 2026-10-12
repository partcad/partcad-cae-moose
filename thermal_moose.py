#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae thermal`, by MOOSE: steady heat conduction through one part.

`temperature:` holds the faces its ports cover at a temperature; `heat:` puts a
power in through them, spread over their area; `convection:` cools them by the
surrounding fluid against `ambient_temperature`. Every other surface is
insulated. MOOSE's heat transfer module solves it, in Celsius -- conduction and
convection do not care where the zero is.
"""

import os
import sys
import tempfile

# `runpy.run_path()` leaves the script's directory off `sys.path`, so the
# sibling module below is not importable without saying where it is.
sys.path.append(os.path.dirname(__file__))
import moose_common as moose  # noqa: E402
from moose_common import WARNING, support  # noqa: E402


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    held, held_records, held_empty = support.condition_groups(request, "temperature", "temperature")
    heated, heated_records, heated_empty = support.condition_groups(request, "heat", "heat")
    cooled, cooled_records, cooled_empty = support.condition_groups(request, "convection", "convection")
    findings.extend(moose.nothing_covered(record) for record in held_empty + heated_empty + cooled_empty)
    if not held_records and not cooled_records:
        raise support.ToolFailed(
            "'thermal:' gives the heat nowhere to go: hold a surface at a 'temperature:', or cool one by "
            "'convection:', on a port that covers a surface of the part. Without one there is no steady state"
        )
    order = int(request.get("mesh_order", 2))

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_MOOSE_CASE") or work
        mesh = moose.mesh_part(request, {**held, **heated, **cooled}, work, order)
        moose.run(
            work,
            _input(request, order, zip(held, held_records), zip(heated, heated_records), zip(cooled, cooled_records)),
        )
        extremes = moose.postprocessors(work)
        nodes = moose.vector(work, "nodes")

    temperature = [value[0] if value else None for value in moose.at_mesh(mesh, nodes, "T")]
    moose.draw(path, mesh, temperature)

    hottest, coldest = extremes["hottest"], extremes["coldest"]
    limit = request.get("max_temperature")
    if limit is not None and hottest > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part reaches %.4g C, past the %.4g C it may." % (hottest, float(limit)),
                "where": "hottest point",
            }
        )
    put_in = sum(float(record["heat"]) for record in heated_records)
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Temperature from %.4g C to %.4g C, %.4g W put in, over %d elements"
            % (coldest, hottest, put_in, mesh.element_count)
        ],
    }


def _input(request, order, held, heated, cooled):
    ambient = float(request.get("ambient_temperature", 20.0))
    conditions = []
    for name, record in held:
        conditions.append(
            "  [hold_%s]\n    type = DirichletBC\n    variable = T\n    value = %.10g\n    boundary = '%s'\n  []"
            % (name, float(record["temperature"]), name)
        )
    for name, record in heated:
        # A power over the area it goes in through: a flux, into the part.
        flux = float(record["heat"]) / (float(record["area"]) * moose.M2_PER_MM2)
        conditions.append(
            "  [heat_%s]\n    type = NeumannBC\n    variable = T\n    value = %.10g\n    boundary = '%s'\n  []"
            % (name, flux, name)
        )
    for name, record in cooled:
        conditions.append(
            "  [cool_%s]\n    type = ConvectiveFluxFunction\n    variable = T\n    coefficient = %.10g\n"
            "    T_infinity = %.10g\n    boundary = '%s'\n  []" % (name, float(record["convection"]), ambient, name)
        )
    return """
[Mesh]
  [file]
    type = FileMeshGenerator
    file = part.msh
  []
[]
[Variables]
  [T]
    order = %(family)s
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
[]
[BCs]
%(conditions)s
[]
[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = %(conductivity).10g
  []
[]
[Postprocessors]
  [hottest]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
  [coldest]
    type = NodalExtremeValue
    variable = T
    value_type = min
  []
[]
[VectorPostprocessors]
  [nodes]
    type = NodalValueSampler
    variable = 'T'
    sort_by = id
  []
[]
[Executioner]
  type = Steady
%(direct)s
[]
[Outputs]
  csv = true
[]
""" % {
        "family": "SECOND" if order == 2 else "FIRST",
        "conditions": "\n".join(conditions),
        "conductivity": float(request.get("conductivity", 50.0)),
        "direct": moose.DIRECT,
    }


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
