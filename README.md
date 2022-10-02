# Relate Python

## Motivation

[MyersGroup's Relate](https://github.com/MyersGroup/relate) is written in C++. Frankly, I am not good at C++. So, I decide to write a wrapper over this project. Maybe I will refactor them all in Python/Cython is the performance loss is acceptable. Maybe I will integrate output file handling (e.g. plotting) in this project.

## Install

```console
$ pip install relatepy  # or
$ pip install git+https://github.com/tcztzy/relate_py
```

## Build from source

Build requirements:

1. C++ compiler support C++17
2. Python interpreter and `<Python.h>` headers
3. zlib
4. Any build tools support [PEP 517](https://peps.python.org/pep-0517/), such as [`build`](https://github.com/pypa/build).

CMake is not required since we have Python's build tool chain.

```console
$ git clone https://github.com/tcztzy/relate_py
$ cd relate_py
$ git submodule init
$ # C++20 had a <version> header, so you should remove the version file in gzstream directory
$ rm relate/include/src/gzstream/version
$ pip install build  # Any tools support PEP 517 would work
$ pip -m build
```

## Environments

Set `RELATEPY_RESOURCE_USAGE` to enable the resource usage report (also need `--verbosity=DEBUG`).
