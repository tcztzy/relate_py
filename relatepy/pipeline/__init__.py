import logging
import struct
from pathlib import Path

import click_log

from ..data import RelateData
from ._build_topology import build_topology
from ._combine_sections import combine_sections
from ._finalize import finalize
from ._find_equivalent_branches import find_equivalent_branches
from ._get_branch_length import get_branch_length
from ._paint import paint as paint_pipeline
from .chunk import chunk as chunk_pipeline

__all__ = ("all_pipeline", "chunk_pipeline", "paint_pipeline")
logger = logging.getLogger(__package__)
click_log.basic_config(logger)


def all_pipeline(
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
    if chunk_index is not None:
        logger.info(f"  chunk {chunk_index}")
        fmt = "ii"
        N, L = struct.unpack(
            fmt,
            (output / f"parameters_c{chunk_index}.bin").read_bytes()[
                : struct.calcsize(fmt)
            ],
        )
        start_chunk = end_chunk = chunk_index
        memory_size = None
    else:
        logger.info(
            "\n  ".join(
                (
                    "Using:",
                    str(haps),
                    str(sample),
                    str(genetic_map),
                    f"with mu = {mutation_rate} and "
                    + (
                        f"2Ne = {effective_population_size}"
                        if coal is None
                        else f"coal = {coal}"
                    ),
                )
            )
        )
        chunk_pipeline(
            haps, sample, genetic_map, output, dist, use_transitions, memory_limit
        )
        fmt = "iiid"
        N, L, end_chunk, memory_size = struct.unpack(
            fmt, (output / "parameters.bin").read_bytes()[: struct.calcsize(fmt)]
        )
        end_chunk -= 1
        start_chunk = 0
    data = RelateData(number_of_sequences=N, number_of_alleles=L)
    logger.info(
        f"Read {data.number_of_sequences} haplotypes with "
        f"{data.number_of_alleles} SNPs per haplotype."
    )
    if memory_size is not None:
        logger.info(f"Expected minimum memory usage: {memory_size}Gb.")

    for c in range(start_chunk, end_chunk + 1):
        logger.info(f"Starting chunk {c} of {end_chunk}.")
        fmt = "i"
        (num_sections,) = struct.unpack_from(
            fmt,
            (output / f"parameters_c{c}.bin").read_bytes()[: struct.calcsize("iii")],
            struct.calcsize("ii"),
        )
        num_sections -= 1
        paint_pipeline(output=output, chunk_index=c, theta=theta, rho=rho)
        build_topology(
            output=output,
            chunk_index=c,
            first_section=0,
            last_section=num_sections - 1,
            effective_population_size=effective_population_size,
            theta=theta,
            rho=rho,
            seed=seed,
            ancestral_state=ancestral_state,
            sample_ages=sample_ages,
        )
        find_equivalent_branches(output=output, chunk_index=c)
        get_branch_length(
            output=output,
            chunk_index=c,
            first_section=0,
            last_section=num_sections - 1,
            mutation_rate=mutation_rate,
            effective_population_size=effective_population_size,
            coal=coal,
            seed=seed,
            sample_ages=sample_ages,
        )
        combine_sections(
            output=output,
            chunk_index=c,
            effective_population_size=effective_population_size,
        )
    if chunk_index is None:
        finalize(output=output, sample_ages=sample_ages, annotation=annotation)
    logger.info("Done.")
