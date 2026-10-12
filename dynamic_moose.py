#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae dynamic`, by MOOSE: the natural frequencies of one part, held as it is.

`fix:` holds the faces its ports cover still; the rest is left to ring. The
generalized eigenproblem K u = omega^2 M u is solved by SLEPc through MOOSE's
`Eigenvalue` executioner, shift-and-invert about zero, which finds the lowest
`modes` frequencies -- the way MOOSE's own modal analysis test poses it, with
the mass as a reaction term tagged for the eigen side.
"""

import math
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
    fixed, fixed_records, fixed_empty = support.condition_groups(request, "fix", "fix")
    findings.extend(moose.nothing_covered(record) for record in fixed_empty)
    if not fixed_records:
        raise support.ToolFailed(
            "'dynamic:' holds nothing, so the part is free and its lowest modes are rigid-body motion at 0 Hz: "
            "name what holds it under 'fix:', on a port that covers a surface of the part"
        )
    modes = max(1, int(request.get("modes", 6)))
    order = int(request.get("mesh_order", 2))

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_MOOSE_CASE") or work
        mesh = moose.mesh_part(request, fixed, work, order)
        moose.run(work, _input(request, fixed, modes, order))
        eigenvalues = moose.vector(work, "omega_squared")
        nodes = moose.vector(work, "nodes")

    squares = sorted(row["eigen_values_real"] for row in eigenvalues)[:modes]
    frequencies = [math.sqrt(max(value, 0.0)) / (2.0 * math.pi) for value in squares]
    if not frequencies:
        raise support.ToolFailed("MOOSE found no modes")
    moved = [value[0] if value else None for value in moose.at_mesh(mesh, nodes, "moved")]
    peak = max([value for value in moved if value is not None], default=0.0) or 1.0
    moose.draw(path, mesh, [None if value is None else value / peak for value in moved], low=0.0)

    limit = request.get("min_frequency")
    if limit is not None and frequencies[0] < float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part's first natural frequency is %.4g Hz, under the %.4g Hz it has to stay above: "
                "it will resonate with anything driving it near that." % (frequencies[0], float(limit)),
                "where": "first natural frequency",
            }
        )
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Natural frequencies %s Hz, over %d elements"
            % (", ".join("%.5g" % value for value in frequencies), mesh.element_count)
        ],
    }


def _input(request, fixed, modes, order):
    family = "SECOND" if order == 2 else "FIRST"
    density = float(request.get("density", 7850))
    blocks = []
    for axis in "xyz":
        blocks.append(
            "  [mass_%s]\n    type = ADCoefReaction\n    variable = disp_%s\n    extra_vector_tags = 'eigen'\n"
            "    coefficient = %.10g\n  []" % (axis, axis, -density)
        )
        blocks.append(
            "  [stiffness_%s]\n    type = StressDivergenceTensors\n    variable = disp_%s\n    component = %d\n  []"
            % (axis, axis, "xyz".index(axis))
        )
    conditions = []
    for axis in "xyz":
        conditions.append(
            "  [hold_%s]\n    type = DirichletBC\n    variable = disp_%s\n    value = 0\n    boundary = %s\n  []"
            % (axis, axis, moose.names(fixed))
        )
        conditions.append(
            "  [hold_%s_eigen]\n    type = EigenDirichletBC\n    variable = disp_%s\n    boundary = %s\n  []"
            % (axis, axis, moose.names(fixed))
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
[Kernels]
%(kernels)s
[]
[BCs]
%(conditions)s
[]
[AuxVariables]
  [moved]
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
  [strain]
    type = ComputeSmallStrain
  []
[]
[Executioner]
  type = Eigenvalue
  solve_type = KRYLOVSCHUR
  which_eigen_pairs = SMALLEST_MAGNITUDE
  n_eigen_pairs = %(modes)d
  n_basis_vectors = %(basis)d
  petsc_options_iname = '-st_type -eps_target -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = 'sinvert 0 lu mumps'
  eigen_tol = 1e-9
[]
[Problem]
  type = EigenProblem
  active_eigen_index = 0
[]
[VectorPostprocessors]
  [omega_squared]
    type = Eigenvalues
    execute_on = TIMESTEP_END
  []
  [nodes]
    type = NodalValueSampler
    variable = 'moved'
    sort_by = id
    execute_on = TIMESTEP_END
  []
[]
[Outputs]
  [csv]
    type = CSV
    file_base = case_out
    execute_on = TIMESTEP_END
  []
[]
""" % {
        "family": family,
        "kernels": "\n".join(blocks),
        "conditions": "\n".join(conditions),
        "youngs": float(request.get("youngs_modulus", 2.1e11)),
        "poisson": float(request.get("poissons_ratio", 0.30)),
        "modes": modes,
        "basis": max(2 * modes + 1, modes + 4),
    }


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
