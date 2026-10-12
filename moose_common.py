#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""Finding MOOSE, running an input file, and reading what it wrote.

Everything the four analysis scripts share that is about MOOSE rather than about
the physics.

**Which MOOSE.** An application built with every physics module -- the one INL
publishes as the `moose` conda package and in its `idaholab/moose` image, where
it answers to `combined-opt`, `moose-opt` and `moose`. A MOOSE *application* of
one's own works too, as long as it was built with the solid mechanics, heat
transfer and Navier-Stokes modules this package's input files use.

**Where it is.** The executable `PARTCAD_MOOSE` names, then each of those three
names on `PATH`, then in the sandbox's own conda environment, then in the usual
places.

**What comes back.** Not MOOSE's own Exodus output, which needs a netCDF reader
to open: every input file here asks for its numbers as CSV -- the extremes as
postprocessors, the fields at every node as a vector postprocessor sampling
them -- which is plain text, and the node positions in it are what the fields
are joined to gmsh's mesh by.
"""

import csv
import glob
import os
import platform
import shutil
import subprocess
import sys

# The meta-wrapper's directory is on 'sys.path' when PartCAD runs a script, and
# 'cae_support' lives there: what every implementation does around its solver.
import cae_support as support  # noqa: E402

EXECUTABLE_ENV = "PARTCAD_MOOSE"
NAMES = ("combined-opt", "moose-opt", "moose")
DIRECTORIES = ("/opt/moose/bin", "/usr/local/bin", os.path.expanduser("~/miniforge3/envs/moose/bin"))

ERROR = "error"
WARNING = "warning"
INFO = "info"
M2_PER_MM2 = 1e-6
MM_PER_M = 1000.0


def find_moose():
    """The MOOSE application to run input files with."""
    named = os.environ.get(EXECUTABLE_ENV)
    if named:
        if os.path.isfile(named) and os.access(named, os.X_OK):
            return named
        raise support.ToolMissing("%s points at something that is not an executable: %s" % (EXECUTABLE_ENV, named))
    for name in NAMES:
        found = shutil.which(name)
        if found:
            return found
    places = [os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "bin")] + list(DIRECTORIES)
    for directory in places:
        for name in NAMES:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
    raise support.ToolMissing(
        "MOOSE is not available on this %s-%s machine. Two ways to have it: run with a container runtime "
        "available, which gets the image named by 'dockerImage:' in partcad.yaml and carries MOOSE -- or install "
        "it here: 'conda create -n moose -c https://conda.software.inl.gov/public moose' on x86_64 Linux and on "
        "macOS, which provides 'combined-opt' with every physics module; then put it on PATH, or point %s at it. "
        "Searched: %s on PATH, then %s"
        % (platform.system(), platform.machine(), EXECUTABLE_ENV, ", ".join(NAMES), ", ".join(places))
    )


def mesh_part(request, groups, work, order, rest=None):
    """Mesh the part as gmsh ASCII, which libMesh reads with its groups' names."""
    msh = os.path.join(work, "part.msh")
    return support.mesh(
        request, groups, size=float(request.get("mesh_size", 0.05)), order=order, rest=rest, files=[msh]
    )


def run(work, text, threads=None):
    """Write an input file into 'work', run MOOSE on it, and return what it printed.

    MOOSE prints an error as a framed block beginning '*** ERROR ***', and that
    block is the whole of what a user needs to read.
    """
    moose = find_moose()
    path = os.path.join(work, "case.i")
    with open(path, "w") as f:
        f.write(text)
    arguments = [moose, "-i", path, "--color", "off"]
    if threads:
        arguments += ["--n-threads=%d" % int(threads)]
    completed = subprocess.run(arguments, cwd=work, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    log = completed.stdout + completed.stderr
    with open(os.path.join(work, "log.moose"), "w") as f:
        f.write(log)
    if completed.returncode != 0 or "*** ERROR ***" in log:
        lines = log.splitlines()
        start = next((k for k, line in enumerate(lines) if "*** ERROR ***" in line), max(len(lines) - 25, 0))
        raise support.ToolFailed(
            "MOOSE failed (exit %d):\n%s" % (completed.returncode, "\n".join(lines[start : start + 25]))
        )
    return log


def postprocessors(work):
    """The postprocessors' values at the last step: name -> number."""
    path = os.path.join(work, "case_out.csv")
    if not os.path.isfile(path):
        raise support.ToolFailed("MOOSE wrote no postprocessor values: the solve did not finish")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise support.ToolFailed("MOOSE wrote an empty postprocessor table")
    return {key: float(value) for key, value in rows[-1].items()}


def vector(work, name):
    """A vector postprocessor's last table, as a list of rows of numbers by column."""
    found = sorted(glob.glob(os.path.join(work, "case_out_%s_*.csv" % name)))
    if not found:
        raise support.ToolFailed("MOOSE wrote no '%s' table" % name)
    with open(found[-1]) as f:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(f)]


def at_mesh(mesh, rows, *columns):
    """The sampled columns at the gmsh mesh's nodes, joined by position."""
    points = [(row["x"], row["y"], row["z"]) for row in rows]
    at = support.match_points(mesh.coordinates, points, mesh.extent * 1e-6)
    return [tuple(rows[k][column] for column in columns) if k is not None else None for k in at]


def draw(path, mesh, values, low=None):
    """The part's surface coloured by one value per gmsh node, as glTF."""
    index = {int(tag): k for k, tag in enumerate(mesh.node_tags)}
    cells = []
    for _type, _tags, nodes in mesh.volumes:
        cells.extend([index[int(tag)] for tag in element[:4]] for element in nodes)
    numbers = [0.0 if value is None else float(value) for value in values]
    support.write_glb(path, mesh.coordinates, support.outer_triangles(cells, 4), numbers, low=low)


def nothing_covered(record):
    """The finding for a condition whose port covers no surface of the part."""
    return {
        "severity": WARNING,
        "message": (
            "The port covers no surface of the part, so the condition on it did nothing. Place the port on "
            "the face it means, or draw its interface with a sketch that lies on one."
        ),
        "where": record.get("port") or record.get("interface_label") or "a port",
    }


def names(groups):
    """A MOOSE list of boundary names."""
    return "'%s'" % " ".join(groups)


DIRECT = """
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
"""
