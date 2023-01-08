import pathlib
from setuptools import setup, Extension
from Cython.Build import cythonize

cwd = pathlib.Path("")
relate = cwd / "relate" / "include"
src = relate / "src"
gzstream = src / "gzstream"
pipeline = relate / "pipeline"

relate_sources = [
    str(src / f"{cpp}.cpp")
    for cpp in (
        "filesystem plot fast_log collapsed_matrix data sample fast_painting "
        "anc mutations tree_builder branch_length_estimator anc_builder "
        "tree_comparer"
    ).split()
]

extensions = [
    Extension(
        "relatepy.data",
        ["relatepy/data.py", str(gzstream / "gzstream.cpp")] + relate_sources,
        include_dirs=[str(dir) for dir in [src, gzstream]],
        libraries=["z"],
        language="c++",
        extra_compile_args=["-std=c++17"],  # use C++ 17 for cpp_locals
        extra_link_args=["-std=c++17"],
    ),
    *[
        Extension(
            f"relatepy.pipeline._{p}",
            [f"relatepy/pipeline/_{p}.py"] + relate_sources,
            include_dirs=[str(dir) for dir in [src, gzstream, pipeline]],
            libraries=["z"],
            language="c++",
            extra_compile_args=["-std=c++17"],  # use C++ 17 for cpp_locals
            extra_link_args=["-std=c++17"],
        ) for p in ("build_topology", "combine_sections", "finalize", "find_equivalent_branches", "get_branch_length")
    ],
]

setup(
    packages=["relatepy"],
    package_data={
        "relatepy": [
            "../relate/include/**/*.cpp",
            "../relate/include/**/*.hpp",
            "options.cpp",
        ]
    },
    ext_modules=cythonize(extensions),
    zip_safe=False,
)
