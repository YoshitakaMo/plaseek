#!/usr/bin/env python3
# %%
from pathlib import Path
from absl import logging
import pandas as pd
import subprocess
import shutil
from absl import app
from absl import flags
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.set_verbosity(logging.INFO)

flags.DEFINE_string(
    "parallel_binary_path",
    shutil.which("parallel"),
    "Path to the parallel executable.",
)
flags.DEFINE_string(
    "tblastn_binary_path", shutil.which("tblastn"), "Path to the tblastn executable."
)
flags.DEFINE_string(
    "foldseek_binary_path", shutil.which("foldseek"), "Path to the Foldseek executable."
)
flags.DEFINE_string(
    "db_path",
    None,
    "Path to the target DNA sequence database.",
)
flags.DEFINE_string(
    "foldseek_tsvfile",
    None,
    "Path to the foldseek m8-formatted file.",
)
flags.DEFINE_string(
    "outfile_path",
    None,
    "Path to the output file.",
)

FLAGS = flags.FLAGS


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
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    retcode = process.wait()
    if retcode != 0:
        raise RuntimeError(
            f"tblastn failed with return code {retcode}.\n" f"stderr:\n{stderr.decode()}"
        )
    with open(outfile, "w") as fh:
        fh.write(stdout.decode())


def collect_pident_plasmid(infile: str, outfile: Path, pident_threshold: float = 98.0) -> None:
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


def main(argv):
    if len(argv) > 1:
        raise app.UsageError("Too many command-line arguments.")
    foldseek_tsvfile = FLAGS.foldseek_tsvfile
    foldseekhits = parse_tsvfile(foldseek_tsvfile)
    logging.info("creating foldseekhits_nodup.fasta")
    with open("foldseekhit_nodups.fasta", "w") as fh:
        SeqIO.write(remove_duplicates(foldseekhits), fh, "fasta")

    run_tblastn(
        tblastn_binary_path=FLAGS.tblastn_binary_path,
        parallel_binary_path=FLAGS.parallel_binary_path,
        db=FLAGS.db_path,
        input_fasta="foldseekhit_nodups.fasta",
        outfile="intermediate.tsv",
    )
    collect_pident_plasmid("intermediate.tsv", FLAGS.outfile_path, pident_threshold=98.0)


if __name__ == "__main__":
    app.run(main)
