from pathlib import Path
import os
import subprocess
from typing import Union
from requests import get, post
from time import sleep
import logging
import tarfile

logger = logging.getLogger(__name__)


def run_foldseek_locally(
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
    pdbfile: str | Path,
    output_directory: str | Path,
    mode="3diaa",
    database: list = ["afdb50", "mgnify_esm30"],
) -> tuple[Path, list[str]]:
    """Retrieve Foldseek results from the webserver.
    Args:
        pdbfile: Path to input pdb file.
        m8file: Path to output m8 file.
        mode: Mode of Foldseek. "3diaa" or "tmalign".
        database: Database of Foldseek.
    Returns:

    """
    # reference : https://github.com/soedinglab/MMseqs2-App/blob/master/docs/api_example.py
    if mode not in ["3diaa", "tmalign"]:
        raise ValueError(f"mode should be 3diaa or tmalign, but got {mode}.")
    if not Path(pdbfile).exists():
        raise FileNotFoundError(f"{pdbfile} not found.")
    file_name = Path(pdbfile).stem
    output_directory = Path(output_directory)
    result_file = output_directory.joinpath(f"result_{file_name}.tar.gz")
    if result_file.exists():
        logger.info(f"{result_file} already exists. Skip running Foldseek.")
    else:
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
            status = get(
                "https://search.foldseek.com/api/ticket/" + ticket["id"]
            ).json()
            logger.info(f"current status is {status['status']}")
            if status["status"] == "ERROR":
                raise RuntimeError(f"Foldseek failed with error: {status['error']}")

            # wait a short time between poll requests
            sleep(5)
            # if status['status'] is not "COMPLETE" repeat = True
            repeat = status["status"] != "COMPLETE"

        download = get(
            "https://search.foldseek.com/api/result/download/" + ticket["id"],
            stream=True,
        )

        with open(result_file, "wb") as fd:
            for chunk in download.iter_content(chunk_size=128):
                fd.write(chunk)

    return result_file, database


def write_merged_m8file(resulttargzfile: Path, mergedm8file: Path) -> Path:
    """Merge m8 files in the {result_file}"""
    with open(mergedm8file, "w") as outfile:
        with tarfile.open(resulttargzfile, "r") as tarf:
            logger.info(f"Extracting m8 files from {resulttargzfile}.")
            for member in tarf.getmembers():
                if member.name.endswith(".m8"):
                    logger.info(f"Found {member.name}")
                    f = tarf.extractfile(member)
                    if f is not None:
                        outfile.write(f.read().decode())
    return Path(mergedm8file)


def run_foldseek_webserver(
    input: Path,
    output_directory: Path,
    mode="3diaa",
    database=["afdb50", "mgnify_esm30"],
) -> Path:
    """run Foldseek webserver and return the path to the merged m8 file.
    Args:
        input: Path to input pdb file.
        mode: Mode of Foldseek. "3diaa" or "tmalign".
        database: Database of Foldseek.
    """
    assert Path(input).exists(), f"{input} not found."
    assert input.suffix == ".pdb", f"input file should be .pdb, but got {input.suffix}."
    resulttargzfile, database = retrieve_foldseek_results(
        input, output_directory, mode, database
    )
    mergedm8file = write_merged_m8file(
        resulttargzfile, Path(output_directory.joinpath(f"{input.stem}.m8"))
    )
    return mergedm8file
