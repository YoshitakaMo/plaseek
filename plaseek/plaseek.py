#!/usr/bin/env python3
# %%
from pathlib import Path
from absl import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import pandas as pd
import subprocess
import shutil
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Union

logging.set_verbosity(logging.INFO)


def check_binaries_available(
    parallel_binary_path: str, tblastn_binary_path: str, foldseek_binary_path: str
) -> None:
    """Check if binaries are available."""
    if not Path(parallel_binary_path).exists():
        raise FileNotFoundError("foldseek not found. Please set PATH to foldseek.")
    if not Path(tblastn_binary_path).exists():
        raise FileNotFoundError("tblastn not found. Please set PATH to tblastn.")
    if not Path(foldseek_binary_path).exists():
        raise FileNotFoundError("parallel not found. Please set PATH to parallel.")


def run_foldseek(
    pdbfile: Union[str, Path],
    foldseek_binary_path: Union[str, Path],
    foldseek_db_path: str,
    outtsvfile: Union[str, Path],
    dbtype: str = "afdb50",
    outfmt: str = "query,target,qstart,pident,fident,nident,qend,qlen,tstart,tend,tlen,alnlen,evalue,bits,qseq,tseq,qaln,taln,qcov,tcov,taxid,taxname,taxlineage,qtmscore,ttmscore,alntmscore,prob",
) -> None:
    """Run Foldseek.
    Args:
        pdbfile: Path to input pdb file.
        foldseek_binary_path: Path to Foldseek binary.
        foldseek_db_path: Path to Foldseek database.
    Raises:
        FileNotFoundError: If binary_path is not found.
    """
    if not Path(foldseek_binary_path).exists():
        raise FileNotFoundError(f"{foldseek_binary_path} not found.")
    if not Path(foldseek_db_path).exists():
        raise FileNotFoundError(f"{foldseek_db_path} not found.")
    cmd = (
        f"{foldseek_binary_path} easy-search {pdbfile} {foldseek_db_path}/{dbtype} {outtsvfile} "
        f"tmp --alignment-type 2 --db-load-mode 2 --max-seqs 1000 -e 10 -s 9.5 --threads 4 "
        f"--prefilter-mode 1 --cluster-search 1 --tmscore-threshold 0.3 --format-mode 4 --remove-tmp-files 1 "
        f"--format-output '{outfmt}'"
    )
    logging.info(f"Launching subprocess: {cmd}")
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    retcode = process.wait()
    if retcode != 0:
        raise RuntimeError(
            f"Foldseek failed with return code {retcode}.\n"
            f"stderr:\n{stderr.decode()}"
        )


def parse_tsvfile(file: Union[str, Path]) -> list[SeqRecord]:
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
    outfmt: str = "6 qaccver saccver pident length qstart qend sstart send qseq sseq qframe sframe evalue bitscore",
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
    if not Path(f"{db}.ndb").exists():
        raise FileNotFoundError(f"{db} not found.")
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


def filtering_by_pident(
    infile: str,
    pident_threshold: float = 98.0,
    sort_values: str = "saccver",
) -> pd.DataFrame:
    """Collect plasmid accession ID from tblastn output file.
    The pident value should be greater than or equal to 98.0 (default).
    """
    df = pd.read_csv(infile, sep="\t", header=None)
    # qaccver saccver pident length qstart qend sstart send qseq sseq qframe sframe evalue bitscore
    df.columns = [
        "qaccver",
        "saccver",
        "pident",
        "length",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qseq",
        "sseq",
        "qframe",
        "sframe",
        "evalue",
        "bitscore",
    ]
    df_filtered = df[df["pident"] >= pident_threshold]
    # sort by sort_values (default: saccver)
    df_filtered_sorted = df_filtered.sort_values(by=[f"{sort_values}"])
    return df_filtered_sorted


def get_unique_seq_regions_from_df(df: pd.DataFrame) -> pd.DataFrame:
    """Get unique sequence regions from a dataframe."""
    df_unique = df.drop_duplicates(subset=["saccver", "sstart", "send"])
    return df_unique


def run_blastdbcmd(
    blastdbcmd_binary_path: str | None,
    parallel_binary_path: str | None,
    db: str,
    df: pd.DataFrame,
) -> None:
    """
    Run blastdbcmd to extract the sequence of the hit region.
    Use GNU parallel to accelerate the process.
    args:
        blastdbcmd_binary_path: Path to the blastdbcmd binary.
        parallel_binary_path: Path to the parallel binary.
        db: Path to the database file.
        df: Dataframe containing the hit regions.
    For example:
    parallel blastdbcmd -db ../db/P476DB -entry {df["sccver"]} -range {df["sstart"]}-{df["send]}
    """
    if blastdbcmd_binary_path is None:
        blastdbcmd_binary_path = str(shutil.which("blastdbcmd"))
    if parallel_binary_path is None:
        parallel_binary_path = str(shutil.which("parallel"))

    with tempfile.NamedTemporaryFile(delete=False, mode="w") as fh:
        df_unique = get_unique_seq_regions_from_df(df)
        # write df_unique["saccver"], df_unique["sstart"], and df_unique["send"] to a temporary file.
        df_unique.to_csv(
            fh,
            sep="\t",
            columns=["saccver", "sstart", "send"],
            header=False,
            index=False,
        )
        tmpfile = fh.name
    cmd = f"parallel -a {tmpfile} --colsep '\\t' 'blastdbcmd -db {db} -entry {{1}} -range {{2}}-{{3}}' ::: {{}}"
    logging.info(f"Launching subprocess: {cmd}")
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    retcode = process.wait()
    if retcode != 0:
        raise RuntimeError(
            f"blastdbcmd failed with return code {retcode}.\n"
            f"stderr:\n{stderr.decode()}"
        )


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
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        help="Path to the input file. pdb or foldseek tsv file are acceptable.",
    )
    parser.add_argument(
        "--foldseek-db-path",
        type=str,
        default=os.getenv("FOLDSEEKDB"),
        help="Path to foldseek database.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 0.0.1",
    )
    tblastn_group = parser.add_argument_group("tblastn arguments", "")
    tblastn_group.add_argument(
        "--target-sequence-db-path",
        type=str,
        default=None,
        help="Path to the target sequence database file for tblastn.",
    )
    tblastn_group.add_argument(
        "--evalue",
        type=float,
        default=1e-100,
        help="E-value threshold for tblastn.",
    )
    tblastn_group.add_argument(
        "--block",
        type=int,
        default=1500,
        help="Block size for tblastn.",
    )
    parser.add_argument(
        "-o",
        "--outfile-path",
        type=str,
        default=None,
        help="Path to the output file.",
    )

    args = parser.parse_args()

    check_binaries_available(
        args.parallel_binary_path, args.tblastn_binary_path, args.foldseek_binary_path
    )

    input = Path(args.input)
    # TODO: available for mmCIF format.
    if input.suffix == ".pdb":
        foldseek_tsvfile = Path(f"{input.stem}.tsv")
        run_foldseek(
            pdbfile=input,
            foldseek_binary_path=args.foldseek_binary_path,
            foldseek_db_path=args.foldseek_db_path,
            outtsvfile=foldseek_tsvfile,
        )
    elif input.suffix == ".tsv":
        foldseek_tsvfile = input
    else:
        raise ValueError("Invalid input file suffix: the suffix must be .pdb or .tsv.")

    foldseekhits = parse_tsvfile(foldseek_tsvfile)

    nodup_fasta = f"{input.stem}_nodup.fasta"
    with open(nodup_fasta, "w") as fh:
        SeqIO.write(remove_duplicates(foldseekhits), fh, "fasta")

    tblastn_result = f"{input.stem}_tblastn.tsv"
    run_tblastn(
        tblastn_binary_path=args.tblastn_binary_path,
        parallel_binary_path=args.parallel_binary_path,
        db=args.target_sequence_db_path,
        input_fasta=nodup_fasta,
        outfile=tblastn_result,
        evalue=args.evalue,
        block=args.block,
    )

    # hit_df = filtering_by_pident(tblastn_result, pident_threshold=98.0)


if __name__ == "__main__":
    main()
