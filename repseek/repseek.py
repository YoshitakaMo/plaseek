#!/usr/bin/env python3
# %%
from pathlib import Path
from absl import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.set_verbosity(logging.INFO)


def parse_tsvfile(file: str) -> list[SeqRecord]:
    """Parse Foldseek result file in TSV format."""
    tsvfile = Path(file)
    df = pd.read_csv(tsvfile, sep="\t", header=0)
    foldseekhits = [
        SeqRecord(
            id=row["target"],
            name=row["target"],
            seq=Seq(row["tseq"]),
            description=f"pident={row['pident']}, taxid={row['taxid']}, taxname={row['taxname']}, "
            f"taxlineage={row['taxlineage']}",
        )
        for _, row in df.iterrows()
    ]
    return foldseekhits


def remove_duplicates(hits: list[SeqRecord]) -> list[SeqRecord]:
    """Remove duplicate sequences from a fasta file.

    Args:
        hits (list[SeqRecord])]):
    Returns:
        list[SeqRecord]: List of SeqRecord objects.
    """
    seen = set()
    nodups = []
    for hit in hits:
        if hit.seq not in seen:
            nodups.append(hit)
            seen.add(hit.seq)
    return nodups


def run_tblastn(
    tblastn_binary_path: str,
    parallel_binary_path: str,
    input_fasta: str,
    db: str,
    outfile: str,
    evalue: float = 1e-100,
    block: int = 2000,
    outfmt: str = "6 qaccver saccver pident length evalue bitscore",
    max_target_seqs: int = 10000,
) -> None:
    """Python wrapper for tblastn
    Args:
        tblastn_binary_path: Path to tblastn binary.
        input_fasta_path: Path to input fasta file.
        db: Path to database file.
        evalue: E-value threshold.
        block: Block size.
        outfmt: Output format.
        max_target_seqs: Maximum number of target sequences.
    Raises:
        FileNotFoundError: If binary_path or db is not found.
    """
    # Requires GNU parallel.
    # https://www.biostars.org/p/76009/
    if not Path(tblastn_binary_path).exists():
        raise FileNotFoundError(f"{tblastn_binary_path} not found.")
    if not Path(parallel_binary_path).exists():
        raise FileNotFoundError(f"{parallel_binary_path} not found.")
    cmd = (
        f"cat {input_fasta} | {parallel_binary_path} --block {block} --recstart '>' --pipe {tblastn_binary_path} "
        f"-evalue {evalue} -db {db} -max_target_seqs {max_target_seqs} -outfmt \\'{outfmt}\\' -query -"
    )
    logging.info(f"Launching subprocess: {cmd}")
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    retcode = process.wait()
    if retcode != 0:
        raise RuntimeError(
            f"tblastn failed with return code {retcode}.\n"
            f"stderr:\n{stderr.decode()}"
        )
    with open(outfile, "w") as fh:
        fh.write(stdout.decode())


def collect_pident_plasmid(
    infile: str, outfile: Path, pident_threshold: float = 98.0
) -> None:
    """Collect plasmid accession ID from tblastn output file.
    The pident value should be greater than or equal to 98.0 (default).
    """
    df = pd.read_csv(infile, sep="\t", header=None)
    df.columns = [
        "qaccver",
        "saccver",
        "pident",
        "length",
        "evalue",
        "bitscore",
    ]
    df_filtered = df[df["pident"] >= pident_threshold]
    with open(outfile, "w") as fh:
        fh.write("\n".join(df_filtered["saccver"].unique()))


def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    binary_group = parser.add_argument_group("binary arguments", "")
    binary_group.add_argument(
        "--parallel_binary_path",
        type=str,
        default=shutil.which("parallel"),
        help="Path to the parallel executable.",
    )
    binary_group.add_argument(
        "--tblastn_binary_path",
        type=str,
        default=shutil.which("tblastn"),
        help="Path to the tblastn executable.",
    )
    binary_group.add_argument(
        "--foldseek_binary_path",
        type=str,
        default=shutil.which("foldseek"),
        help="Path to the Foldseek executable.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )
    parser.add_argument(
        "--db_path",
        type=str,
        default=None,
        help="Path to the database file.",
    )
    parser.add_argument(
        "--foldseek_tsvfile",
        type=str,
        default=None,
        help="Path to the foldseek m8-formatted file.",
    )
    parser.add_argument(
        "--outfile_path",
        type=str,
        default=None,
        help="Path to the output file.",
    )

    args = parser.parse_args()
    foldseek_tsvfile = args.foldseek_tsvfile
    foldseekhits = parse_tsvfile(foldseek_tsvfile)
    logging.info("creating foldseekhits_nodup.fasta")
    with open("foldseekhit_nodups.fasta", "w") as fh:
        SeqIO.write(remove_duplicates(foldseekhits), fh, "fasta")

    run_tblastn(
        tblastn_binary_path=args.tblastn_binary_path,
        parallel_binary_path=args.parallel_binary_path,
        db=args.db_path,
        input_fasta="foldseekhit_nodups.fasta",
        outfile="intermediate.tsv",
    )
    collect_pident_plasmid("intermediate.tsv", args.outfile_path, pident_threshold=98.0)


if __name__ == "__main__":
    main()
