from pathlib import Path
from absl import logging
import subprocess
import shutil
from typing import Optional


def run_tblastn(
    input_fasta: str | Path,
    db: str,
    outfile: str | Path,
    block: int = 3000,
    tblastn_binary_path: Optional[str] = shutil.which("tblastn"),
    parallel_binary_path: Optional[str] = shutil.which("parallel"),
    evalue: float = 1e-50,
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
    if tblastn_binary_path is None:
        raise FileNotFoundError(f"{tblastn_binary_path} not found.")
    if parallel_binary_path is None:
        raise FileNotFoundError(f"{parallel_binary_path} not found.")
    if not Path(f"{db}.ndb").exists():
        raise FileNotFoundError(f"{db}.ndb not found.")

    input_fasta = str(input_fasta)
    db = str(db)
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

    header = outfmt[2:].replace(" ", "\t")
    with open(outfile, "w") as fh:
        fh.write(f"{header}\n")
        fh.write(stdout.decode())
