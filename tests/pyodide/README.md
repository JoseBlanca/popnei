# popnei under pyodide

pyodide is CPython built for WebAssembly, which runs in a browser tab and
under node. popnei is installed in it as a wheel whose platform tag is
`pyemscripten`, built from the same core crate and the same binding crate
as the native wheel, and what this directory checks is that the wheel
builds, installs and answers with the version of the core.

From the root of the repository, from a clean checkout:

    scripts/build_pyodide_wheel.sh
    (cd tests/pyodide && npm install)
    node tests/pyodide/smoke.mjs

The first command prints the path of the wheel it left in `dist/`; the
third prints `popnei.__version__` as pyodide gives it and exits with an
error when that is not the version in `[workspace.package]` of the
`Cargo.toml` of the repository. Neither `dist/` nor `node_modules/` is in
git.

## What has to be installed

Outside the repository, on the machine that builds. "What has to be in
place" of `docs/plans/vcf-to-blocks.md` has the commands that installed
them here.

- The Rust target `wasm32-unknown-emscripten`, `rustup target add
  wasm32-unknown-emscripten`. cargo names it when it is not installed.
- emsdk, the installer of emscripten, in `~/devel/emsdk`, with 5.0.3
  activated in it. emscripten is the C compiler and linker that compiles to
  WebAssembly, and the CPython of pyodide 314.0.7 was compiled with that
  version of it. The build script sources `emsdk_env.sh`, which puts `emcc`
  on the PATH together with the node and the Python that emsdk downloaded
  for it: `emcc` refuses to run under a Python older than 3.10, and the
  only `python3` of the PATH on this Mac is Apple's 3.9. Set `EMSDK` to the
  directory of another emsdk to build with that one.
- pyodide-build in the virtual environment `~/devel/pyodide-venv`, with the
  cross build environment of pyodide 314.0.7 installed in it, which is
  where the headers and the libraries of the CPython of pyodide are. Set
  `PYODIDE_VENV` to the directory of another virtual environment.
- node, for the smoke test, which takes pyodide itself from npm.

The build script stops with a message when emsdk or pyodide-build is not
where it looks, when the host Python of pyodide-build is the free threaded
build of CPython, and when the version of `emcc` is not the one pyodide
asks for.

## Moving to another version of pyodide

Three things name the version, and the three have to agree, because a wheel
built against one build of the CPython of pyodide does not load in another:
the cross build environment, `pyodide xbuildenv install <version>`; the
emscripten that was activated in `~/devel/emsdk`; and the `pyodide`
dependency of `package.json` beside this file, which is the runtime the
smoke test loads.

The build script reads the version of emscripten from `pyodide config get`
and checks `emcc` against it, so changing the cross build environment and
the emsdk is enough on that side; the npm dependency is the one line to
change by hand.

## Two traps

Both were found when the trial crate in Rust that `docs/rust_core.md`
reports was built for pyodide.

- The host Python that pyodide-build runs under must not be the free
  threaded build of CPython, the one without the global interpreter lock,
  which `uv venv --python 3.14` picks on this Mac. With it pyodide-build
  writes the sysconfig data of emscripten under `python3.14t`, the `t` of
  which marks that build, where nothing looks for it. The build script
  stops when the host Python of pyodide-build is that one, and
  `.python-version`, the file uv reads to choose the interpreter of the
  `.venv` of popnei, names 3.14.5 and not 3.14 for the same reason.
- `pyodide config get rust_toolchain` answers 1.93.0, older than the
  `rust-version = "1.98"` of the workspace, and that is not a problem:
  pyodide-build installs and switches to that toolchain in
  `pyodide build-recipes`, which builds the packages of the pyodide
  distribution from its own recipes, and not in `pyodide build`, which is
  what this script runs. The wheel here is compiled by whichever `cargo`
  is on the PATH, 1.98.0 on this machine.

## What runs where

The smoke test writes the wheel into the file system that emscripten gives
the code running inside pyodide, which is not the one of node, and installs
it from there with micropip, the package installer of pyodide. The `emfs:`
prefix of the path is what asks micropip for a file of that file system
instead of a package of an index. micropip itself is not in the npm package
of pyodide: the first run downloads it from the CDN of pyodide and caches
it under `node_modules/`, so that first run needs the network and the later
ones do not.
