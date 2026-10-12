#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae cfd`, by MOOSE: steady, laminar, incompressible flow through one part.

The part is read as the **fluid volume**. `load:` names an inlet, held at the
pressure that is its force over the area its port covers; `outlet:` names where
the flow leaves, at the reference pressure; every other surface is a no-slip
wall. Where an inlet or an outlet faces along an axis, the flow through it is
also held normal to it, which is what a developed flow through a cut does.

MOOSE's Navier-Stokes module, the finite element formulation with PSPG and SUPG
stabilization on linear tetrahedra, solved by Newton with a direct factorization
-- the input MOOSE's own pressure-driven channel test uses, in three dimensions.
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

TURBULENT_REYNOLDS = 2300.0


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    density = float(request.get("density", 1.204))
    viscosity = float(request.get("viscosity", 1.82e-5))
    inlets, inlet_records, inlet_empty = support.condition_groups(request, "load", "inlet")
    outlets, outlet_records, outlet_empty = support.condition_groups(request, "outlet", "outlet")
    findings.extend(moose.nothing_covered(record) for record in inlet_empty + outlet_empty)
    if not inlet_records:
        raise support.ToolFailed("'cfd:' drives nothing: name an inlet under 'load:', on a port that covers a surface")
    if not outlet_records:
        raise support.ToolFailed(
            "'cfd:' names no outlet, so the flow has nowhere to go and the problem has no downstream pressure "
            "reference: name where it leaves under 'outlet:', on a port that covers a surface"
        )
    pressures = [float(r["load"]) / (float(r["area"]) * moose.M2_PER_MM2) for r in inlet_records]

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_MOOSE_CASE") or work
        mesh = moose.mesh_part(request, {**inlets, **outlets}, work, 1, rest="walls")
        moose.run(work, _input(request, inlets, inlet_records, pressures, outlets, outlet_records, density, viscosity))
        nodes = moose.vector(work, "nodes")

    sampled = moose.at_mesh(mesh, nodes, "vel_x", "vel_y", "vel_z", "p")
    flow = sum(
        support.flow_through(
            mesh, name, lambda k: sampled[k][:3] if sampled[k] else None, support.outward_normal(record)
        )
        for name, record in zip(outlets, outlet_records)
    )
    moose.draw(path, mesh, [value[3] if value else None for value in sampled])

    outlet_area = sum(float(record["area"]) for record in outlet_records) * moose.M2_PER_MM2
    mean_speed = flow / outlet_area if outlet_area else 0.0
    peak_speed = max(math.sqrt(row["vel_x"] ** 2 + row["vel_y"] ** 2 + row["vel_z"] ** 2) for row in nodes)
    reynolds = density * abs(mean_speed) * math.sqrt(4.0 * outlet_area / math.pi) / viscosity

    if reynolds > TURBULENT_REYNOLDS:
        findings.append(
            {
                "severity": WARNING,
                "message": "The flow is likely turbulent (Reynolds number about %.0f, past %.0f), and this analysis "
                "is laminar." % (reynolds, TURBULENT_REYNOLDS),
            }
        )
    limit = request.get("max_velocity")
    if limit is not None and peak_speed > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The flow reaches %.3g m/s, past the %.3g m/s it is allowed." % (peak_speed, float(limit)),
                "where": "peak speed",
            }
        )
    minimum = request.get("min_flow_rate")
    if minimum is not None and flow < float(minimum):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part passes %.4g m^3/s, short of the %.4g m^3/s it has to." % (flow, float(minimum)),
                "where": "flow rate",
            }
        )
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Flow rate %.4g m^3/s (%.4g L/min), mean outlet speed %.4g m/s, peak speed %.4g m/s, Reynolds number "
            "%.3g, driven by %s, over %d elements"
            % (
                flow,
                flow * 60000.0,
                mean_speed,
                peak_speed,
                reynolds,
                ", ".join("%.4g Pa" % value for value in pressures),
                mesh.element_count,
            )
        ],
    }


def _normal_only(record):
    """The velocity functions that hold a flow normal to an axis-aligned face, or None."""
    normal = support.outward_normal(record)
    if normal is None:
        return None
    axis = max(range(3), key=lambda k: abs(normal[k]))
    if abs(normal[axis]) < 0.999:
        return None
    return " ".join("NA" if k == axis else "0" for k in range(3))


def _input(request, inlets, inlet_records, pressures, outlets, outlet_records, density, viscosity):
    velocity_boundaries, velocity_functions = ["walls"], ["0 0 0"]
    for name, record in list(zip(inlets, inlet_records)) + list(zip(outlets, outlet_records)):
        functions = _normal_only(record)
        if functions:
            velocity_boundaries.append(name)
            velocity_functions.append(functions)
    pressure_boundaries = list(inlets) + list(outlets)
    pressure_functions = ["%.10g" % value for value in pressures] + ["0"] * len(outlets)
    return """
[Mesh]
  [file]
    type = FileMeshGenerator
    file = part.msh
  []
[]
[Modules]
  [IncompressibleNavierStokes]
    equation_type = steady-state
    velocity_boundary = '%(vb)s'
    velocity_function = '%(vf)s'
    pressure_boundary = '%(pb)s'
    pressure_function = '%(pf)s'
    density_name = rho
    dynamic_viscosity_name = mu
    integrate_p_by_parts = false
    order = FIRST
    pspg = true
    supg = true
  []
[]
[Materials]
  [fluid]
    type = GenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '%(rho).10g %(mu).10g'
  []
[]
[VectorPostprocessors]
  [nodes]
    type = NodalValueSampler
    variable = 'vel_x vel_y vel_z p'
    sort_by = id
  []
[]
[Executioner]
  type = Steady
%(direct)s
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-14
  nl_max_its = %(iterations)d
[]
[Outputs]
  csv = true
[]
""" % {
        "vb": " ".join(velocity_boundaries),
        "vf": "  ".join(velocity_functions),
        "pb": " ".join(pressure_boundaries),
        "pf": " ".join(pressure_functions),
        "rho": density,
        "mu": viscosity,
        "direct": moose.DIRECT,
        "iterations": int(request.get("iterations", 30)),
    }


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
