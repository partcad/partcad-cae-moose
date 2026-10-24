# partcad-cae-moose

Static stress, modal, steady heat transfer and laminar flow analysis for
[PartCAD](https://github.com/partcad/partcad), by [MOOSE](https://mooseframework.inl.gov/), INL's multiphysics
framework. It implements all four of PartCAD's analyses:

```shell
pc cae fea -i //pub/feature/cae/moose:fea :bracket
pc cae dynamic -i //pub/feature/cae/moose:dynamic :bracket
pc cae thermal -i //pub/feature/cae/moose:thermal :heatsink
pc cae cfd -i //pub/feature/cae/moose:cfd :duct
```

or name it for one of a part's scenarios with `implementation:` in it, or for everything with the
`cae<Analysis>Implementation` user configuration options.

## What each analysis reads

Each analysis reads its own section of the part, a map of named scenarios: `pc cae <analysis>` runs each, and
`-f <scenario>` one. Every key names interfaces the part implements, and is applied to the **faces each port
covers**: PartCAD imprints a port's sketch on the part before this package runs; a port with no sketch covers
the face it lies on.

| analysis | section keys | MOOSE |
| --- | --- | --- |
| `fea` | `fix:`, `load:`, `gravity:` | the solid mechanics module, small strain; `DirichletBC` on every displacement; `NeumannBC`, the force over the faces' area, along the direction the part states (`[x, y, z]` in its frame) or along -Z for a magnitude; `Gravity` for the part's own weight, under the `gravity:` it states or 9.81 m/s² along -Z |
| `dynamic` | `fix:` | the same stiffness, the mass tagged for the eigen side, `Eigenvalue` with SLEPc shift-and-invert for the lowest `modes` |
| `thermal` | `temperature:`, `heat:`, `convection:` | `HeatConduction`; `DirichletBC` (C); `NeumannBC`, the power over the faces' area; `ConvectiveFluxFunction` against `ambient_temperature`. Every other surface is insulated |
| `cfd` | `load:`, `outlet:`, `fix:` | the Navier-Stokes module, finite elements with PSPG and SUPG; an inlet held at its force over its area, an outlet at the reference pressure, every other surface a no-slip wall |

### What the numbers mean

A `load` value is a force, and PartCAD turns it into newtons before this package sees it. Write a number and a
unit -- `n`, `mn`, `kn` or `newton` for a force, `mg`, `g`, `kg`, `ton`, `tonne`, `lb` or `pound` for a mass,
matched case-insensitively, with or without a space and a plural `s` -- or a **bare number, which is a mass in
kilograms**, weighed at 9.8 N/kg: `hook: 5` is 49 N. A direction is three of those, each read the same way, so
`[0, 0, -5]` is the same 49 N straight down, and `[0 N, -2.7 N, -155.2 N]` is newtons as written. What a
simulation measured is in newtons, so write the unit: `[0, -2.7, -155.2]` would be read as 1.5 kN.
`cfd` reads its `load:` the same way: the newtons that push the fluid in, spread over the inlet's area.

`gravity:` is an acceleration, m/s² along each of the part's axes -- a number, or a string saying so
(`9.81 m/s^2`). The part's own weight is its mass, `density` times its volume, under that acceleration, in
place of the 9.81 m/s² along -Z that `self_weight` applies otherwise; it is never a `load`.

The model is `<part>.<analysis>.<scenario>.glb` -- coloured by von Mises stress, by the first mode's
displacement, by temperature, or by pressure. A measurement line goes to the log, and **findings** fail
`pc test`: a stress past the yield strength (or past it divided by `safety_factor`), a displacement past
`max_displacement`, a first natural frequency under `min_frequency`, a temperature past `max_temperature`, a
flow fast enough to be turbulent or past `max_velocity` or short of `min_flow_rate`, and any condition whose
port covers no surface.

### The pipeline

gmsh meshes the part into tetrahedra -- quadratic, with quadratic unknowns, for the structural and thermal
analyses; linear for the flow -- each condition's faces a named group (PartCAD's `cae_support.mesh()`), and
libMesh reads the groups by name. The script writes the input file, and MOOSE writes its answers as CSV: the
extremes as postprocessors, and every node's position and values as a vector postprocessor, which is joined to
gmsh's nodes by position for the model. No Exodus reader is needed.

### Measured

`examples/feature_cae` in the PartCAD repository has a cantilever and a pipe whose answers are formulas.
Through `pc cae`, with MOOSE 2026.09.25 from INL's conda channel:

| analysis | closed form | this package |
| --- | --- | --- |
| `fea`: tip deflection, 100 N on a 100x10x10 mm steel bar | 0.2016 mm (with shear) | 0.2005 mm |
| `dynamic`: first bending frequency | 815.4 Hz (Euler-Bernoulli) | 813.25 Hz, twice -- the square section's two bending planes |
| `thermal`: hottest point, 1 W at the tip, root at 20 C | 40.0 C | 40.0 C |
| `cfd`: flow rate, 5 mN through a 10x100 mm bore of oil | 1.5625 mL/s | 1.544 mL/s |

Elmer and Code_Aster, on the same mesh, give the same displacement, frequencies and temperatures to every
printed digit.

## What it needs

**A container runtime on x86_64, or MOOSE installed.**

With a container runtime PartCAD runs the analyses in the image `partcad.yaml` names: PartCAD's Python base
image with INL's `moose` conda package -- `combined-opt`, every physics module in one application -- in a prefix
of its own, and Debian's gmsh. **amd64 only**, and large: INL publishes MOOSE for 64-bit x86 Linux and macOS
only, and it brings PETSc, libMesh and VTK with it.

Without one, PartCAD runs them in a conda or venv sandbox on the host, and MOOSE has to be there. It is found,
in this order, where `PARTCAD_MOOSE` points, as `combined-opt`, `moose-opt` or `moose` on `PATH`, in the
sandbox's own conda environment, and in `/opt/moose/bin`, `/usr/local/bin`, `~/miniforge3/envs/moose/bin`. The
easiest way to have it:

```shell
conda create -n moose -c https://conda.software.inl.gov/public moose
```

which works on x86_64 Linux and on macOS (Intel and Apple Silicon). A MOOSE application of your own works too,
as long as it was built with the solid mechanics, heat transfer and Navier-Stokes modules.

The daemon is what runs the analysis, so a variable set in your shell reaches it only when the daemon is
started from that shell: `pc daemon stop` first.

### The pin

`dockerImage:` names the image by a hash of what it is built from -- the Dockerfile, `verify_image.py` and the
Python version. Run `./image-tag.sh` after changing any of them and update every pin; the `Image` workflow
refuses to build if they disagree, and publishes `<tag>-amd64` on a push to `main`.

## Tests

```shell
pytest
```

beside an installed `partcad`. The input files, the CSV readers and the messages a machine without MOOSE gets
are checked everywhere; the cantilever is analysed three ways and compared with its closed forms where MOOSE,
gmsh and build123d are all present.

To keep a case for inspection, set `PARTCAD_MOOSE_CASE` to an empty directory before the daemon starts.

## License

Apache License 2.0. MOOSE is LGPL-2.1 and is not part of this package: the image installs INL's published
conda build of it.
