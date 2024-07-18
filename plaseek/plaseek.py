#!/usr/bin/env python3
# %%
from pathlib import Path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import pandas as pd
import shutil
import time
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Union
import plaseek.tools.foldseek
import plaseek.tools.blastdbcmd
import plaseek.tools.tblastn
from plaseek.tools.utils import setup_logging
import plaseek.tools.utils

logger = logging.getLogger(__name__)


def check_binaries_available(
    parallel_binary_path: str, tblastn_binary_path: str
) -> None:
    """Check if binaries are available."""
    if not Path(parallel_binary_path).exists():
        raise FileNotFoundError("parallel not found. Please set PATH to GNU parallel.")
    if not Path(tblastn_binary_path).exists():
        raise FileNotFoundError("tblastn not found. Please set PATH to tblastn.")


def filtering_m8file(
    file: Union[str, Path], eval_threshold: float = 1e-10
) -> pd.DataFrame:
    """Filtering Foldseek result file in M8 format.
    filter by e-value < 1e-10
    The header is "query,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,prob,evalue,bits,qlen,tlen,qaln,taln,tca,tseq,taxid,taxname".
    """
    m8file = Path(file)
    df = pd.read_csv(
        m8file,
        sep="\t",
        names=[
            "query",
            "theader",
            "pident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "prob",
            "evalue",
            "bits",
            "qlen",
            "tlen",
            "qaln",
            "taln",
            "tca",
            "tseq",
            "taxid",
            "taxname",
        ],
    )
    df_filtered = df[df["evalue"] < eval_threshold]
    return df_filtered


def remove_duplicates(hits: pd.DataFrame) -> list[SeqRecord]:
    """Remove duplicate sequences from a fasta file.

    Args:
        hits (pd.DataFrame):
    Returns:
        list[SeqRecord]: List of SeqRecord objects.
    """
    seen = set()
    nodups = []
    for _, row in hits.iterrows():
        if row["tseq"] not in seen:
            nodups.append(
                SeqRecord(
                    id=row["theader"],
                    name=row["theader"],
                    seq=Seq(row["tseq"]),
                    description=f"pident={row['pident']}, evalue={row['evalue']} taxid={row['taxid']}, taxname={row['taxname']}, ",
                )
            )
            seen.add(row["tseq"])
    return nodups


def filtering_by_pident(
    infile: str | Path,
    minpident_threshold: float = 98.0,
    maxpident_threshold: float = 100.0,
    sort_values: str = "saccver",
) -> pd.DataFrame:
    """Collect plasmid accession ID from tblastn output file.
    The pident value should be greater than or equal to 98.0 (default).
    """
    df = pd.read_csv(infile, sep="\t", header=0)

    df_filtered = df[
        (df["pident"] >= minpident_threshold) & (df["pident"] <= maxpident_threshold)
    ]
    # sort by sort_values (default: saccver)
    df_filtered_sorted = df_filtered.sort_values(by=[f"{sort_values}"])
    return df_filtered_sorted


def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    binary_group = parser.add_argument_group("binary arguments", "")
    binary_group.add_argument(
        "--parallel-binary-path",
        type=str,
        default=shutil.which("parallel"),
        help="Path to the parallel executable.",
    )
    binary_group.add_argument(
        "--tblastn-binary-path",
        type=str,
        default=shutil.which("tblastn"),
        help="Path to the tblastn executable.",
    )
    binary_group.add_argument(
        "--foldseek-binary-path",
        type=str,
        default=shutil.which("foldseek"),
        help="Path to the Foldseek executable.",
    )
    binary_group.add_argument(
        "--blastdbcmd-binary-path",
        type=str,
        default=shutil.which("blastdbcmd"),
        help="Path to the blastdbcmd executable.",
    )
    webserver_group = parser.add_argument_group("Foldseek webserver arguments", "")
    webserver_group.add_argument(
        "--alignmode",
        type=str,
        choices=["3diaa", "tmalign"],
        default="3diaa",
        help="Mode for Foldseek. 3diaa or tmalign.",
    )
    webserver_group.add_argument(
        "--foldseek_web_database",
        type=list[str],
        choices=[
            "pdb100",
            "afdb50",
            "afdb-swissprot",
            "afdb-proteome",
            "cath50",
            "mgnify_esm30",
            "gmgcl_id",
        ],
        default=["afdb50", "mgnify_esm30"],
        help="Sequence database for Foldseek webserver.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        required=True,
        help="Path to the input file. pdb or foldseek tsv file are acceptable.",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        default=None,
        help="Path to the output directory.",
    )
    parser.add_argument(
        "-f",
        "--foldseek-db-path",
        type=str,
        default=os.getenv("FOLDSEEKDB"),
        help="Path to foldseek database.",
    )
    parser.add_argument(
        "-t",
        "--target-sequence-db-path",
        type=str,
        default=None,
        required=True,
        help="Path to the target sequence database file for tblastn.",
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 0.0.3",
    )
    foldseek_group = parser.add_argument_group("Foldseek arguments", "")
    foldseek_group.add_argument(
        "--foldseek_evalue_threshold",
        type=float,
        default=1e-20,
        help="evalue threshold for Foldseek.",
    )
    tblastn_group = parser.add_argument_group("tblastn arguments", "")
    tblastn_group.add_argument(
        "--tblastn_evalue_threshold",
        type=float,
        default=1e-50,
        help="tblastn evalue threshold.",
    )
    tblastn_group.add_argument(
        "--tblastn_minpident_threshold",
        type=float,
        default=50.0,
        help="min pident threshold for tblastn. use 0.0-100.0.",
    )
    tblastn_group.add_argument(
        "--tblastn_maxpident_threshold",
        type=float,
        default=100.0,
        help="max pident threshold for tblastn. use 0.0-100.0.",
    )
    tblastn_group.add_argument(
        "--block",
        type=int,
        default=3000,
        help="Block size for tblastn.",
    )
    tblastn_group.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=4,
        help="Number of parallel jobs.",
    )

    args = parser.parse_args()
    check_binaries_available(args.parallel_binary_path, args.tblastn_binary_path)

    input = Path(args.input)

    output_directory = args.output_directory
    if output_directory is None:
        output_directory = Path(os.getcwd()).joinpath(f"{input.stem}")
    else:
        output_directory = Path(output_directory)
    os.makedirs(output_directory, exist_ok=True)
    setup_logging(Path(output_directory).joinpath("log.txt"))

    foldseek_m8file = output_directory.joinpath(f"{input.stem}.m8")
    if input.suffix == ".pdb":
        if not Path(foldseek_m8file).exists():
            if args.foldseek_binary_path is None or args.foldseek_db_path is None:
                # use Foldseek webserver
                logger.info(
                    f"foldseek_binary_path is not provided. Use Foldseek webserver for {input}."
                )
                foldseek_m8file = plaseek.tools.foldseek.run_foldseek_webserver(
                    input, output_directory
                )
            elif args.foldseek_binary_path is None:
                raise ValueError(
                    "Please set PATH to Foldseek by --foldseek-binary-path."
                )
            else:
                if args.foldseek_db_path is None:
                    raise ValueError(
                        "Please set PATH to Foldseek database by --foldseek-db-path."
                    )
                foldseek_m8file = Path(f"{input.stem}.m8")
                plaseek.tools.foldseek.run_foldseek_locally(
                    pdbfile=input,
                    foldseek_binary_path=args.foldseek_binary_path,
                    foldseek_db_path=args.foldseek_db_path,
                    outtsvfile=foldseek_m8file,
                    jobs=args.jobs,
                )
        else:
            logger.info(
                f"{foldseek_m8file} already exists. Skip using Foldseek webserver."
            )
    elif input.suffix == ".m8":
        foldseek_m8file = input
    else:
        raise ValueError("Invalid input file suffix: the suffix must be .pdb or .m8.")

    filtered_foldseekhits: pd.DataFrame = filtering_m8file(
        foldseek_m8file, eval_threshold=args.foldseek_evalue_threshold
    )

    nodup_fasta = Path(output_directory.joinpath(f"{input.stem}_nodup.fasta"))
    with open(nodup_fasta, "w") as fh:
        SeqIO.write(remove_duplicates(filtered_foldseekhits), fh, "fasta")

    tblastn_result = Path(output_directory.joinpath(f"{input.stem}_tblastn.tsv"))
    start = time.perf_counter()
    plaseek.tools.tblastn.run_tblastn(
        db=args.target_sequence_db_path,
        input_fasta=nodup_fasta,
        outfile=tblastn_result,
        block=args.block,
        tblastn_binary_path=args.tblastn_binary_path,
        parallel_binary_path=args.parallel_binary_path,
        evalue=args.tblastn_evalue_threshold,
    )
    duration = time.perf_counter() - start
    logger.info(f"tblastn took {duration:.2f} seconds for {input.stem}")

    filtered_tblastn = filtering_by_pident(
        tblastn_result,
        minpident_threshold=args.tblastn_minpident_threshold,
        maxpident_threshold=args.tblastn_maxpident_threshold,
    )
    # TODO: if no hits, skip blastdbcmd
    if filtered_tblastn.empty:
        logger.info(f"No hits found for {input.stem}. Skip blastdbcmd. No results.")
    else:
        logger.info(f"Found {len(filtered_tblastn)} hits for {input.stem}.")
        start = time.perf_counter()
        plaseek.tools.blastdbcmd.run_blastdbcmd(
            blastdbcmd_binary_path=args.blastdbcmd_binary_path,
            parallel_binary_path=args.parallel_binary_path,
            outfile=Path(output_directory.joinpath(f"{input.stem}_result.fasta")),
            db=args.target_sequence_db_path,
            df=filtered_tblastn,
            jobs=args.jobs,
        )
        duration = time.perf_counter() - start
        logger.info(f"blastdbcmd took {duration:.2f} seconds for {input.stem}")

        plaseek.tools.utils.write_resultfile(
            blastcmdresult=output_directory.joinpath(f"{input.stem}_result.fasta"),
            tblastnresult=tblastn_result,
            outputfile=output_directory.joinpath(f"{input.stem}_output.tsv"),
        )


if __name__ == "__main__":
    main()
