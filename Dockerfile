# The runtime the MOOSE analyses run best in.
#
# It carries what pip cannot: a MOOSE application with every physics module,
# which is a native program, and Debian's gmsh. Everything else the scripts
# import is declared in `partcad.yaml` as `pythonRequirements` and installed by
# PartCAD into the environment it creates -- this image is a *runtime*, and the
# package has to work in a conda or venv sandbox on a host with MOOSE installed
# too.
#
# amd64 only: INL publishes MOOSE for 64-bit x86 Linux and for macOS, and not
# for 64-bit ARM Linux. It is also large -- MOOSE brings PETSc, libMesh, VTK and
# a compiler for its just-in-time features -- which is the price of a solver that
# does all four analyses in one program.
ARG PYTHON_VERSION=3.11
ARG ARCH=amd64
FROM ghcr.io/partcad/partcad-container-python:py${PYTHON_VERSION}-${ARCH}

LABEL partcad.image="1"
LABEL org.opencontainers.image.source="https://github.com/partcad/partcad-cae-moose"

USER root
ENV DEBIAN_FRONTEND=noninteractive

# The mesher, from Debian. gmsh's Python API is one ctypes module wrapping
# `libgmsh.so`, installed for Debian's interpreter; the one in this image is
# python.org's, so the module is copied where both can see it and named on
# PYTHONPATH -- PartCAD runs the scripts in a virtual environment it creates at
# run time, which sees neither site-packages.
RUN apt-get update \
  && apt-get install --yes --no-install-recommends python3-gmsh \
  && rm -rf /var/lib/apt/lists/* \
  && mkdir -p /opt/pc-site \
  && cp /usr/lib/python3/dist-packages/gmsh.py /opt/pc-site/gmsh.py
ENV PYTHONPATH=/opt/pc-site

# MOOSE from INL's conda channel, into a prefix of its own: the package pins
# hundreds of exact versions, its own Python among them, so it cannot share an
# environment with anything else and must not shadow the interpreter PartCAD
# runs in. micromamba is fetched by version, checked against the hash
# conda-forge publishes, and removed with its cache in the same layer.
RUN set -eu; \
    MICROMAMBA_VERSION=2.9.0; \
    MOOSE_VERSION=2026.09.25; \
    SHA256=8761c382127e6363bd9e0a2451aa3ef90d071a79133f736e2f759a3bf13040dd; \
    python3 -c 'import hashlib, io, os, sys, tarfile, urllib.request; url, want = sys.argv[1], sys.argv[2]; blob = urllib.request.urlopen(url).read(); got = hashlib.sha256(blob).hexdigest(); got == want or sys.exit("%s hashes to %s, not the pinned %s" % (url, got, want)); open("/usr/local/bin/micromamba", "wb").write(tarfile.open(fileobj=io.BytesIO(blob), mode="r:bz2").extractfile("bin/micromamba").read()); os.chmod("/usr/local/bin/micromamba", 0o755)' \
      "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MICROMAMBA_VERSION}-0.tar.bz2" "${SHA256}"; \
    micromamba create --yes --root-prefix /tmp/mamba --prefix /opt/moose \
      --channel https://conda.software.inl.gov/public --channel conda-forge "moose=${MOOSE_VERSION}"; \
    rm -rf /tmp/mamba /usr/local/bin/micromamba

# Where `moose_common.find_moose()` looks even when PATH says nothing. A symlink
# rather than an edit to PATH, because PartCAD builds a sandbox inside this
# container and a sandbox rebuilds PATH.
RUN ln -s /opt/moose/bin/combined-opt /usr/local/bin/combined-opt

ENV PC_CONTAINER_ALLOWED_COMMANDS='{"python": "/usr/local/bin/python3", "combined-opt": "/usr/local/bin/combined-opt"}'

# Two checks: the base image's own, and this package's, which solves a small
# heat problem with the solver from *this* image, the way the analyses will.
COPY verify_image.py /tmp/verify_image.py
RUN python3 /pc/verify.py \
  && python3 /tmp/verify_image.py \
  && rm /tmp/verify_image.py

USER pc
