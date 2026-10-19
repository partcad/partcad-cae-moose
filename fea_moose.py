#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae fea`, by MOOSE: a linear static stress analysis of one part.

`fix:` holds the faces its ports cover still in every direction; `load:` spreads
its force evenly over the faces its ports cover: along the direction the part
states, a load written as `[x, y, z]` in its own frame, and along -Z, where a
weight goes, for a load written as a magnitude. Gravity acts on the part's own
mass too: the `gravity:` the part states, in its own frame -- what a simulation
found it under, `g - a` -- or 9.81 m/s^2 along -Z unless `self_weight` is off.
MOOSE's solid mechanics module solves it, small strain, on quadratic
tetrahedra with quadratic displacements.
"""

import os
import sys
import tempfile

# `runpy.run_path()` leaves the script's directory off `sys.path`, so the
# sibling module below is not importable without saying where it is.
sys.path.append(os.path.dirname(__file__))
import moose_common as moose  # noqa: E402
from moose_common import ERROR, WARNING, support  # noqa: E402


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    fixed, fixed_records, fixed_empty = support.condition_groups(request, "fix", "fix")
    loaded, loaded_records, loaded_empty = support.condition_groups(request, "load", "load")
    findings.extend(moose.nothing_covered(record) for record in fixed_empty + loaded_empty)
    if not fixed_records:
        raise support.ToolFailed(
            "'fea:' holds nothing, so the part is free to float: name at least one interface under 'fix:', "
            "on a port that covers a surface of the part"
        )

    order = int(request.get("mesh_order", 2))
    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_MOOSE_CASE") or work
        mesh = moose.mesh_part(request, {**fixed, **loaded}, work, order)
        moose.run(work, _input(request, fixed, loaded, loaded_records, order))
        extremes = moose.postprocessors(work)
        nodes = moose.vector(work, "nodes")

    # At the nodes, projected from the elements, which is where the other
    # implementations report it too: an element's own value is the average over
    # it, and reads lower wherever the stress changes fast -- at a clamp.
    peak_stress = max([row["stress"] for row in nodes], default=extremes["peak_stress"])
    peak_displacement = extremes["peak_displacement"] * moose.MM_PER_M
    stress = [value[0] if value else None for value in moose.at_mesh(mesh, nodes, "stress")]
    moose.draw(path, mesh, stress, low=0.0)

    findings.extend(_verdict(peak_stress, peak_displacement, request))
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Peak von Mises stress %.4g MPa, peak displacement %.4g mm, over %d elements"
            % (peak_stress / 1e6, peak_displacement, mesh.element_count)
        ],
    }


def _input(request, fixed, loaded, loaded_records, order):
    density = float(request.get("density", 7850))
    gravity = _gravity_of(request)
    family = "SECOND" if order == 2 else "FIRST"
    held = "\n".join(
        "  [hold_%s]\n    type = DirichletBC\n    variable = disp_%s\n    value = 0\n    boundary = %s\n  []"
        % (axis, axis, moose.names(fixed))
        for axis in "xyz"
    )
    # A force on the port, spread over the area it covers: a traction, so a
    # finer mesh does not load the part harder -- one condition per direction it
    # has a component along.
    pulled = "\n".join(
        "  [pull_%s_%s]\n    type = NeumannBC\n    variable = disp_%s\n    value = %.10g\n    boundary = '%s'\n  []"
        % (name, axis, axis, value / (float(record["area"]) * moose.M2_PER_MM2), name)
        for name, record in zip(loaded, loaded_records)
        for axis, value in zip("xyz", _force_of(record))
        if value != 0.0
    )
    weight = "\n".join(
        "  [weight_%s]\n    type = Gravity\n    variable = disp_%s\n    value = %.10g\n  []" % (axis, axis, value)
        for axis, value in zip("xyz", gravity)
    )
    return """
[Mesh]
  [file]
    type = FileMeshGenerator
    file = part.msh
  []
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[Variables]
  [disp_x]
    order = %(family)s
  []
  [disp_y]
    order = %(family)s
  []
  [disp_z]
    order = %(family)s
  []
[]
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    generate_output = 'vonmises_stress'
    # Linear within each element rather than one value per element, so the
    # projection to the nodes keeps the peak at a clamp instead of averaging
    # it away.
    material_output_order = FIRST
  []
[]
[Kernels]
%(weight)s
[]
[AuxVariables]
  [moved]
  []
  [stress]
  []
[]
[AuxKernels]
  [moved]
    type = VectorMagnitudeAux
    variable = moved
    x = disp_x
    y = disp_y
    z = disp_z
  []
  [stress]
    # At the end, after the stress it projects has been computed: MOOSE works
    # out nodal auxiliary variables before elemental ones, and the von Mises
    # stress is elemental, so in the same pass this would project zeros.
    type = ProjectionAux
    variable = stress
    v = vonmises_stress
    execute_on = final
  []
[]
[BCs]
%(held)s
%(pulled)s
[]
[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = %(youngs).10g
    poissons_ratio = %(poisson).10g
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = %(density).10g
  []
[]
[Postprocessors]
  [peak_displacement]
    type = NodalExtremeValue
    variable = moved
  []
  [peak_stress]
    type = ElementExtremeValue
    variable = vonmises_stress
  []
[]
[VectorPostprocessors]
  [nodes]
    type = NodalValueSampler
    variable = 'moved stress'
    sort_by = id
    execute_on = final
  []
[]
[Executioner]
  type = Steady
%(direct)s
[]
[Outputs]
  [csv]
    type = CSV
    file_base = case_out
    execute_on = 'timestep_end final'
  []
[]
""" % {
        "family": family,
        "weight": weight,
        "held": held,
        "pulled": pulled,
        "youngs": float(request.get("youngs_modulus", 2.1e11)),
        "poisson": float(request.get("poissons_ratio", 0.30)),
        "density": density,
        "direct": moose.DIRECT,
    }


def _force_of(record):
    """The force on one loaded port, newtons in the part's frame: as stated, or its magnitude straight down."""
    force = record.get("force")
    if isinstance(force, (list, tuple)) and len(force) == 3:
        return [float(component) for component in force]
    return [0.0, 0.0, -float(record["load"])]


def _gravity_of(request):
    """The acceleration of free fall the part is weighed under, m/s^2 in its own frame.

    The part's own `gravity:` where it stated one -- a gravity of nought weighs
    nothing -- and otherwise 9.81 m/s^2 along -Z, unless `self_weight` is off.
    """
    stated = request.get("gravity")
    if isinstance(stated, (list, tuple)) and len(stated) == 3:
        return [float(component) for component in stated]
    return [0.0, 0.0, -9.81 if request.get("self_weight", True) else 0.0]


def _verdict(peak_stress, peak_displacement, request):
    """What the numbers mean, as findings: the thresholds the package configured."""
    findings = []
    yield_strength = float(request.get("yield_strength", 2.5e8))
    safety_factor = float(request.get("safety_factor", 1.5))
    allowable = yield_strength / max(safety_factor, 1e-9)
    if peak_stress > yield_strength:
        findings.append(
            {
                "severity": ERROR,
                "message": "The part yields under this load: peak stress %.3g MPa against a yield strength of %.3g MPa."
                % (peak_stress / 1e6, yield_strength / 1e6),
                "where": "peak von Mises stress",
            }
        )
    elif peak_stress > allowable:
        findings.append(
            {
                "severity": WARNING,
                "message": "The part does not reach the safety factor of %.2f: peak stress %.3g MPa against an "
                "allowable %.3g MPa." % (safety_factor, peak_stress / 1e6, allowable / 1e6),
                "where": "peak von Mises stress",
            }
        )
    limit = request.get("max_displacement")
    if limit is not None and peak_displacement > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part moves %.3g mm under this load, past the %.3g mm it is allowed."
                % (peak_displacement, float(limit)),
                "where": "peak displacement",
            }
        )
    return findings


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
