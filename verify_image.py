#!/usr/bin/env python3
#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""Whether this image can actually run the MOOSE analyses.

Run at build time, so that an image which cannot is never published. The base
image's own `verify.py` checks the contract PartCAD talks to a container
through; this checks what this image adds: gmsh imported by *this* interpreter,
and `combined-opt` found where the analyses look for it and solving a bar held
at 0 at one end and 10 at the other, read back the way the analyses read their
answers -- as CSV.
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile

BAR = """
[Mesh]
  [bar]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 1
    nz = 1
  []
[]
[Variables]
  [T]
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
[]
[BCs]
  [cold]
    type = DirichletBC
    variable = T
    boundary = left
    value = 0
  []
  [hot]
    type = DirichletBC
    variable = T
    boundary = right
    value = 10
  []
[]
[Materials]
  [k]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 1
  []
[]
[Postprocessors]
  [hottest]
    type = NodalExtremeValue
    variable = T
  []
[]
[Executioner]
  type = Steady
[]
[Outputs]
  csv = true
[]
"""


def main() -> int:
    failures = []
    try:
        import gmsh

        gmsh.initialize()
        gmsh.finalize()
        print("gmsh imports from this interpreter")
    except Exception as e:
        failures.append("gmsh could not be imported from this interpreter (%s)" % e)

    moose = shutil.which("combined-opt")
    if not moose:
        failures.append("combined-opt is not on PATH")
    else:
        with tempfile.TemporaryDirectory() as work:
            with open(os.path.join(work, "bar.i"), "w") as f:
                f.write(BAR)
            completed = subprocess.run(
                [moose, "-i", "bar.i", "--color", "off"], cwd=work, capture_output=True, text=True
            )
            table = os.path.join(work, "bar_out.csv")
            if completed.returncode != 0 or not os.path.isfile(table):
                failures.append(
                    "combined-opt could not solve a heat problem: %s" % (completed.stdout + completed.stderr)[-800:]
                )
            else:
                with open(table) as f:
                    hottest = float(list(csv.DictReader(f))[-1]["hottest"])
                print("combined-opt solved a bar held at 0 and 10: hottest point %.3f" % hottest)
                if abs(hottest - 10.0) > 1e-6:
                    failures.append("combined-opt's answer to a bar held at 0 and 10 peaks at %g, not 10" % hottest)

    if failures:
        print("\nThis image cannot run the MOOSE analyses:")
        for failure in failures:
            print("  * %s" % failure)
        return 1
    print("\nMOOSE and the mesher are both usable from this interpreter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
