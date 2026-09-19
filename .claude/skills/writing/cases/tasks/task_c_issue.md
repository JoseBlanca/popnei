# Task C: a GitHub issue

While building the Rust spike as a wasm wheel for pyodide you lost forty
minutes to a build failure. Write the GitHub issue for the popnei
repository that records what happened, so that whoever writes popnei's
wasm build script makes it guard against this. Below are your session
notes. They are all you know: where the notes do not hold a fact, do not
supply it. Write the issue, title and body, as markdown, to the output
path you were given.

## Session notes

```
wasm wheel build, spike crate (pyo3 0.29 + maturin), pyodide-build 0.39.0,
xbuildenv pyodide 314.0.7, emscripten 5.0.3, rust stable 1.98

made the venv for pyodide-build with: uv venv --python 3.14
build failed. error was about the emscripten sysconfigdata module not being
found. didn't keep the exact text.

dug around ~40 min. cause: on this machine uv resolves "3.14" to the FREE
THREADED build (3.14t). with that host python pyodide-build installs the
emscripten sysconfigdata under a python3.14t directory, and the build then
looks for it under python3.14. nothing looks in the 3.14t dir.

fix: make the venv with the normal (GIL) 3.14 explicitly. commands are in
spike/README.md in the pyNei repo. after that: builds in 20 s, wheel 193 KB,
micropip installs it in 0.13 s.

don't know: whether uv picks 3.14t on other machines or only here (maybe
because 3.14t was the only 3.14 installed?). didn't check. don't know if
newer pyodide-build handles it.

idea: build script checks sysconfig.get_config_var("Py_GIL_DISABLED") of the
host python and stops with a message if it's 1.
```
