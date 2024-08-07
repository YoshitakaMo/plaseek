import logging
from logging import StreamHandler, Formatter
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from dataclasses import dataclass  # Python 3.7+


def setup_logging(log_file: Path, mode: str = "w") -> None:
    log_file.parent.mkdir(exist_ok=True, parents=True)
    logger = logging.getLogger(__name__)
    if logger.handlers:
        for handler in logger.handlers:
            handler.close()
            logger.removeHandler(handler)
    logger.setLevel(logging.INFO)
    stream_handler = StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    handler_format = Formatter(
        "%(asctime)s %(filename)s:%(lineno)d - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(filename)s:%(lineno)d - %(levelname)s - %(message)s",
        handlers=[stream_handler, logging.FileHandler(log_file, mode=mode)],
    )


def write_resultfile(
    blastcmdresult: str | Path, tblastnresult: str | Path, outputfile: str | Path
):
    """
    Write the result file.
    ----
    blastcmdresult: str | Path
        Path to the result file of blastcmd.
    tblastnresult: str | Path
        Path to the result file of tblastn.
    outputfile: str | Path
        Path to the output file.
    """

    @dataclass
    class BlastcmdHit:
        genbank_id: str
        qseq1: int
        qseq2: int
        dnaseq: str

    @dataclass
    class TblastnHit:
        saccver: str
        sstart: int
        send: int
        sseq: str

    blastcmdhits = []
    with open(blastcmdresult, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))
        for record in records:
            # record.idを:と-で分離してそれぞれを取得
            genbank_id = str(record.id.split(":")[0])
            qseq1 = int(record.id.split(":")[1].split("-")[0])
            qseq2 = int(record.id.split(":")[1].split("-")[1])
            dnaseq = str(record.seq)
            blastcmdhits.append(BlastcmdHit(genbank_id, qseq1, qseq2, dnaseq))

    tblastnhits = []
    with open(tblastnresult, "r") as f:
        df = pd.read_csv(f, sep="\t")
        for idx, row in df.iterrows():
            saccver = str(row["saccver"])
            sstart = int(row["sstart"])
            send = int(row["send"])
            sseq = str(row["sseq"])
            tblastnhits.append(TblastnHit(saccver, sstart, send, sseq))

    with open(outputfile, "w") as f:
        for blastcmdhit in blastcmdhits:
            for tblastnhit in tblastnhits:
                if blastcmdhit.genbank_id == tblastnhit.saccver:
                    if (
                        blastcmdhit.qseq1 == tblastnhit.sstart
                        and blastcmdhit.qseq2 == tblastnhit.send
                    ) or (
                        blastcmdhit.qseq2 == tblastnhit.sstart
                        and blastcmdhit.qseq1 == tblastnhit.send
                    ):
                        if tblastnhit.sstart > tblastnhit.send:
                            # translate from 3' to 5'
                            translated_seq = (
                                Seq(blastcmdhit.dnaseq).reverse_complement().translate()
                            )
                        else:
                            # translate from 5' to 3'
                            translated_seq = Seq(blastcmdhit.dnaseq).translate()
                        f.write(
                            f"{blastcmdhit.genbank_id}\t{tblastnhit.sstart}\t{tblastnhit.send}\t{blastcmdhit.dnaseq}\t{translated_seq}\n"
                        )
