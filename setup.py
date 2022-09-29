import pathlib
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize

cwd = pathlib.Path("")
relate = cwd / "relate"
src = relate / "include" / "src"

extensions = Extension(
    "relatepy.data",
    ["relatepy/data.py"],
    include_dirs=[str(src.parent / "pipeline"), str(src), str(src / "gzstream")],
    library_dirs=[str(relate / "bin")],
    libraries=["relateStatic", "gzstreamStatic", "z"],
    language="c++",
    extra_compile_args=["-std=c++17"],  # use C++ 17 for cpp_locals
    extra_link_args=["-std=c++17"],
)

setup(
    name="relatepy",
    packages=["relatepy"],
    package_data={"relatepy": ["../relate/include/**/*", "../relate/bin/*", "options.cpp"]},
    ext_modules=cythonize(extensions),
    zip_safe=False,
)
