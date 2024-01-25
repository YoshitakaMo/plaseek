from absl import logging
import os
import pandas as pd
import subprocess
import shutil
import tempfile
from typing import Optional


def swap_sstart_and_send(df: pd.DataFrame) -> pd.DataFrame:
    """Swap sstart and send if a row has larger "sstart" value than "send" value.
    args:
        df: tblastn-output dataframe.
    returns:
        sanitized_df: tblastn-output dataframe.
    """
    sanitized_df = df.copy()
    for idx, row in df.iterrows():
        if row["sstart"] > row["send"]:
            sanitized_df.at[idx, "sstart"] = row["send"]
            sanitized_df.at[idx, "send"] = row["sstart"]
    return sanitized_df


def run_blastdbcmd(
    db: str,
    df: pd.DataFrame,
    outfile: str,
    blastdbcmd_binary_path: Optional[str] = shutil.which("blastdbcmd"),
    parallel_binary_path: Optional[str] = shutil.which("parallel"),
    jobs: int | None = os.cpu_count(),
) -> None:
    """
    Run blastdbcmd to extract the sequence of the hit region.
    Use GNU parallel to accelerate the process.
    args:
        blastdbcmd_binary_path: Path to the blastdbcmd binary.
        parallel_binary_path: Path to the parallel binary.
        db: Path to the database file.
        jobs: Number of parallel jobs.
        df: Dataframe containing the hit regions.
    E.g.:
    parallel -j <NCPUS> -a /tmp/foo.txt blastdbcmd -db ../db/P847DB -entry {df["sccver"]} -range {df["sstart"]}-{df["send]}
    """
    if blastdbcmd_binary_path is None:
        blastdbcmd_binary_path = str(shutil.which("blastdbcmd"))
    if parallel_binary_path is None:
        parallel_binary_path = str(shutil.which("parallel"))
    if jobs is None:
        jobs = 4

    # Get unique sequence regions from a dataframe.
    df_unique = df.drop_duplicates(subset=["saccver", "sstart", "send"])
    # Swap sstart and send if a row has larger "sstart" value than "send" value.
    sanitized_df = swap_sstart_and_send(df_unique)
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as fh:
        sanitized_df.to_csv(fh, sep="\t", header=False, index=False)
        tmpfile = fh.name

    cmd = f"parallel -j {jobs} -a {tmpfile} --colsep '\\t' 'blastdbcmd -db {db} -entry {{2}} -range {{7}}-{{8}}' ::: {{}}"
    logging.info(f"Launching subprocess with {jobs} jobs: {cmd}")
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
    with open(outfile, "w") as fh:
        fh.write(stdout.decode())
