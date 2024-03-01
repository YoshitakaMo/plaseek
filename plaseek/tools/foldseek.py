from pathlib import Path
from logging import StreamHandler, Formatter
import os
import subprocess
from typing import Union
from requests import get, post
from time import sleep
import logging
import tarfile


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = StreamHandler()
stream_handler.setLevel(logging.DEBUG)
handler_format = Formatter(
    "%(asctime)s %(filename)s:%(lineno)d - %(levelname)s - %(message)s"
)
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)


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

    logger.info(f"Foldseek is running for {pdbfile}.")
    cmd = (
        f"{foldseek_binary_path} easy-search {pdbfile} {foldseek_db_path}/{dbtype} {outtsvfile} "
        f"tmp --alignment-type 2 --max-seqs 1000 -e 10 -s 9.5 --threads 4 "
        f"--prefilter-mode 1 --cluster-search 1 --tmscore-threshold 0.3 --remove-tmp-files 1 "
        f"--db-load-mode 2 --format-output '{outfmt}'"
    )
    logger.info(f"Launching subprocess: {cmd}")
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


def retrieve_foldseek_results(
    pdbfile: str | Path, mode="3diaa", database: list = ["afdb50", "mgnify_esm30"]
) -> tuple[Path, list[str]]:
    """Retrieve Foldseek results from the webserver.
    Args:
        pdbfile: Path to input pdb file.
        m8file: Path to output m8 file.
        mode: Mode of Foldseek.
        database: Database of Foldseek.
    Returns:

    """
    # reference : https://github.com/soedinglab/MMseqs2-App/blob/master/docs/api_example.py
    if not Path(pdbfile).exists():
        raise FileNotFoundError(f"{pdbfile} not found.")
    file_name = Path(pdbfile).stem

    files = {"q": open(pdbfile, "rb")}
    data = {
        "mode": mode,  # select the mode you want to use
        "database[]": database,  # select the database you want to use
    }

    ticket = post(
        "https://search.foldseek.com/api/ticket", files=files, data=data
    ).json()

    print(ticket)

    # poll until the job was successful or failed
    repeat = True
    while repeat:
        status = get("https://search.foldseek.com/api/ticket/" + ticket["id"]).json()
        logger.info(f"current status is {status['status']}")
        if status["status"] == "ERROR":
            raise RuntimeError(f"Foldseek failed with error: {status['error']}")

        # wait a short time between poll requests
        sleep(5)
        # if status['status'] is not "COMPLETE" repeat = True
        repeat = status["status"] != "COMPLETE"

    download = get(
        "https://search.foldseek.com/api/result/download/" + ticket["id"], stream=True
    )

    result_file = Path(f"result_{file_name}.tar.gz")
    with open(result_file, "wb") as fd:
        for chunk in download.iter_content(chunk_size=128):
            fd.write(chunk)
    return result_file, database


def write_merged_m8file(
    result_file: Path, mergedm8file: str = "alis_merged.m8"
) -> Path:
    """Merge m8 files in the {result_file}.tar.gz."""
    with open(mergedm8file, "w") as outfile:
        with tarfile.open(result_file, "r") as tarf:
            logger.info(f"Extracting m8 files from {result_file}.tar.gz")
            for member in tarf.getmembers():
                if member.name.endswith(".m8"):
                    logger.info(f"Found {member.name}")
                    f = tarf.extractfile(member)
                    if f is not None:
                        outfile.write(f.read().decode())
    return Path(mergedm8file)
