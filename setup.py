import pathlib
from setuptools import setup, Extension
from Cython.Build import cythonize

cwd = pathlib.Path.cwd()
relate = cwd / "relate"
src = relate / "include" / "src"
print(relate)

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
    ext_modules=cythonize(extensions),
    zip_safe=False,
)
