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
from dataclasses import dataclass  # Python 3.7+
from typing import Union
import plaseek.tools.foldseek
import plaseek.tools.blastdbcmd
import plaseek.tools.tblastn
from plaseek.tools.utils import setup_logging

resultblastcmd = "/Users/YoshitakaM/Desktop/AF-A0A166M635-F1-model_v4/AF-A0A166M635-F1-model_v4_result.fasta"
tblastnresult = "/Users/YoshitakaM/Desktop/AF-A0A166M635-F1-model_v4/AF-A0A166M635-F1-model_v4_tblastn.tsv"


def write_resultfile(blastcmdhits, tblastnhits, outputfile):
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
    with open(resultblastcmd, "r") as f:
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


write_resultfile(resultblastcmd, tblastnresult, "/Users/YoshitakaM/Desktop/result.tsv")
# %%
