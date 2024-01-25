from pathlib import Path
from absl import logging
import os
import subprocess
from typing import Union


def run_foldseek(
    pdbfile: Union[str, Path],
    foldseek_binary_path: Union[str, Path],
    foldseek_db_path: str,
    outtsvfile: Union[str, Path],
    jobs: int | None = os.cpu_count(),
    dbtype: str = "afdb50",
    outfmt: str = "query,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,prob,evalue,bits,qlen,tlen,qaln,taln,tca,tseq,taxid,taxname",
) -> None:
    """Run Foldseek. The output m8 format (outfmt) should be compatbile with the webserver of Foldseek (https://search.foldseek.com/search)
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

    logging.info(f"Foldseek is running for {pdbfile}.")
    cmd = (
        f"{foldseek_binary_path} easy-search {pdbfile} {foldseek_db_path}/{dbtype} {outtsvfile} "
        f"tmp --alignment-type 2 --max-seqs 1000 -e 10 -s 9.5 --threads 4 "
        f"--prefilter-mode 1 --cluster-search 1 --tmscore-threshold 0.3 --remove-tmp-files 1 "
        f"--db-load-mode 2 --format-output '{outfmt}'"
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
