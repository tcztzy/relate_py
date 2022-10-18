import logging
from pathlib import Path

import click
import click_log

from .pipeline import all_pipeline, chunk_pipeline, paint_pipeline


class NotRequiredIf(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if = kwargs.pop("not_required_if")
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (
            kwargs.get("help", "")
            + " NOTE: This argument is mutually exclusive with %s"
            % self.not_required_if
        ).strip()
        super().__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        we_are_present = self.name in opts
        other_present = self.not_required_if in opts

        if other_present:
            if we_are_present:
                raise click.UsageError(
                    "Illegal usage: `%s` is mutually exclusive with `%s`"
                    % (self.name, self.not_required_if)
                )
            else:
                self.prompt = None

        return super().handle_parse_result(ctx, opts, args)


BASE_DIR = Path(__file__).parent.parent
logger = logging.getLogger(__package__)
click_log.basic_config(logger)
PathType = click.Path(exists=True, path_type=Path)
GLOBAL_OPTIONS = dict(
    haps=(
        "--haps",
        dict(
            required=True,
            help="Filename of haps file (Output file format of Shapeit).",
            type=PathType,
        ),
    ),
    sample=(
        "--sample",
        dict(
            required=True,
            help="Filename of sample file (Output file format of Shapeit).",
            type=PathType,
        ),
    ),
    genetic_map=(
        "--genetic-map",
        "--map",
        dict(
            required=True,
            help="Genetic map.",
            type=PathType,
        ),
    ),
    output=(
        "--output",
        "-o",
        dict(
            required=True,
            help="Filename of output without file extension.",
            type=click.Path(path_type=Path),
        ),
    ),
    sample_ages=(
        "--sample-ages",
        "--sample_ages",
        dict(
            help="Filename of file containing sample ages (one per line).",
            type=PathType,
        ),
    ),
    mutation_rate=(
        "--mutation-rate",
        "--mutation_rate",
        "-m",
        dict(
            required=True,
            help="Mutation rate.",
            type=float,
        ),
    ),
    effective_populate_size=(
        "--effective-population-size",
        "--effectiveN",
        "-N",
        dict(
            cls=NotRequiredIf,
            help="Effective population size.",
            type=float,
            not_required_if="coal",
        ),
    ),
    dist=(
        "--dist",
        dict(
            help="Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used.",
            type=PathType,
        ),
    ),
    annotation=(
        "--annotation",
        "--annot",
        dict(
            help="Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used.",
            type=PathType,
        ),
    ),
    coal=(
        "--coal",
        dict(
            help="Filename of file containing coalescent rates. If specified, it will overwrite --effective-population-size.",
            type=PathType,
        ),
    ),
    chunk_index=(
        "--chunk-index",
        "--chunk_index",
        dict(
            help="Index of chunk. (Use when running parts of the algorithm on an individual chunk.)",
            type=int,
        ),
    ),
    use_transitions=(
        "--transversion",
        "use_transitions",
        dict(
            is_flag=True,
            default=True,
            help="Only use transversion for bl estimation.",
        ),
    ),
    memory_limit=(
        "--memory-limit",
        "--memory",
        dict(
            default=5,
            help="Approximate memory allowance in GB for storing distance matrices.",
            type=float,
        ),
    ),
)


def global_options(*options, **option_args):
    def wrap(func):
        for o in options + tuple(option_args.keys()):
            if o not in GLOBAL_OPTIONS:
                raise click.BadArgumentUsage(f"Unknown key `{o}`.")
            *args, kwargs = GLOBAL_OPTIONS.get(o)
            _kwargs = option_args.get(o, ({},))
            if isinstance(_kwargs, dict):
                _args = []
            else:
                *_args, _kwargs = _kwargs
            kwargs.update(_kwargs)
            func = click.option(*args, *_args, **kwargs)(func)
        return click_log.simple_verbosity_option(logger)(func)

    return wrap


relate = click.Group(
    help="""\
\b
Relate
 * Authors: Leo Speidel, Marie Forest, Sinan Shi, Simon Myers.
 * Doc:     https://myersgroup.github.io/relate
"""
)


@relate.command
@global_options(
    "haps",
    "sample",
    "genetic_map",
    "output",
    "mutation_rate",
    "effective_populate_size",
    "sample_ages",
    "coal",
    "dist",
    "annotation",
    "chunk_index",
    "use_transitions",
    "memory_limit",
)
@click.option(
    "--theta",
    default=0.001,
    type=float,
)
@click.option(
    "--rho",
    default=1.0,
    type=float,
)
@click.option(
    "--anc_allele_unknown",
    "ancestral_state",
    help="Specify if ancestral allele is unknown.",
    is_flag=True,
    default=True,
)
@click.option(
    "--seed",
    help="Seed for MCMC in branch lengths estimation.",
    type=int,
)
def all(
    haps: Path,
    sample: Path,
    genetic_map: Path,
    output: Path,
    mutation_rate: float,
    effective_population_size: float | None,
    sample_ages: list[float],
    dist: Path | None = None,
    annotation: Path | None = None,
    coal: Path | None = None,
    chunk_index: int | None = None,
    use_transitions: bool = True,
    memory_limit: float = 5,
    theta: float = 0.001,
    rho: float = 1,
    ancestral_state: bool = True,
    seed: int | None = None,
):
    all_pipeline(
        haps,
        sample,
        genetic_map,
        output,
        mutation_rate,
        effective_population_size,
        sample_ages,
        dist,
        annotation,
        coal,
        chunk_index,
        use_transitions,
        memory_limit,
        theta,
        rho,
        ancestral_state,
        seed,
    )


@relate.command
@global_options(
    "haps", "sample", "genetic_map", "output", "dist", "use_transitions", "memory_limit"
)
def chunk(
    haps: Path,
    sample: Path,
    genetic_map: Path,
    output: Path,
    dist: Path | None = None,
    use_transitions: bool = True,
    memory_limit: float = 5.0,
) -> None:
    """Chunk the input data."""
    try:
        chunk_pipeline(
            haps, sample, genetic_map, output, dist, use_transitions, memory_limit
        )
    except FileExistsError as e:
        raise click.FileError(
            e.filename,
            "Directory already exists. Relate will use this directory to store temporary files.",
        )


@relate.command
@global_options("output", chunk_index=dict(required=True))
@click.option(
    "--theta",
    default=0.001,
    type=float,
)
@click.option(
    "--rho",
    default=1.0,
    type=float,
)
def paint(output: Path, chunk_index: int, theta: float, rho: float):
    paint_pipeline(output, chunk_index, theta, rho)


if __name__ == "__main__":
    relate()
