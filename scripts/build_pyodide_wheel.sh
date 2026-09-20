#!/usr/bin/env bash
# Builds the wheel of popnei for pyodide and leaves it in dist/.
# tests/pyodide/README.md says what has to be installed on the machine, how
# to run this script and how to run the smoke test that installs the wheel.

set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

# Neither emscripten nor pyodide-build is a dependency of the repository:
# they are installed on the machine that builds, where "What has to be in
# place" of docs/plans/vcf-to-blocks.md puts them. Set either variable to
# build with another copy.
emsdk_dir=${EMSDK:-$HOME/devel/emsdk}
pyodide_venv=${PYODIDE_VENV:-$HOME/devel/pyodide-venv}

if [ ! -f "$emsdk_dir/emsdk_env.sh" ]; then
    echo "no emsdk in $emsdk_dir: set EMSDK to the directory it is in" >&2
    exit 1
fi
if [ ! -x "$pyodide_venv/bin/pyodide" ]; then
    echo "no pyodide-build in $pyodide_venv: set PYODIDE_VENV to its venv" >&2
    exit 1
fi

# emcc, and the node and the Python that emsdk downloaded for it: emcc
# refuses to run under a Python older than 3.10, and the only python3 of
# the PATH on this Mac is Apple's 3.9.
export EMSDK_QUIET=1
# shellcheck source=/dev/null
. "$emsdk_dir/emsdk_env.sh"

export PATH="$pyodide_venv/bin:$PATH"

# The free threaded build of CPython as the host Python of pyodide-build
# leaves the sysconfig data of emscripten under python3.14t, where nothing
# looks for it, and the wheel is built without a single error against the
# wrong interpreter. Section 3.2 of docs/rust_core.md has the trap.
if [ "$("$pyodide_venv/bin/python" -c 'import sys; print(sys._is_gil_enabled())')" != "True" ]; then
    echo "the host Python of $pyodide_venv is the free threaded one" >&2
    exit 1
fi

# The cross build environment of pyodide was compiled with one version of
# emscripten and links against no other, so the version is read from it and
# not written here.
wanted_emscripten=$(pyodide config get emscripten_version)
found_emscripten=$(emcc -dumpversion)
if [ "$found_emscripten" != "$wanted_emscripten" ]; then
    echo "emcc is $found_emscripten and pyodide asks for $wanted_emscripten:" \
        "install that one in $emsdk_dir and activate it" >&2
    exit 1
fi

# An older wheel left in dist/ would still be there if this build failed,
# and the smoke test, which takes the one wheel it finds, would install it
# and pass.
rm -f "$repo_root"/dist/popnei-*-pyemscripten*.whl

cd "$repo_root"
pyodide build

shopt -s nullglob
wheels=("$repo_root"/dist/popnei-*-pyemscripten*.whl)
if [ ${#wheels[@]} -ne 1 ]; then
    echo "pyodide build left ${#wheels[@]} wheels with a pyemscripten tag" \
        "in dist/, and it has to leave one" >&2
    exit 1
fi
echo "built ${wheels[0]}"
