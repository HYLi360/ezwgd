"""
Program frontend.
"""

import pandas as pd
import tempfile

import subprocess
import os

def gffread(
        gff3_file_path: str,
        fna_file_path: str,
        cds_out: str,
        pep_out: str,
        abs_gffread_path: str = 'gffread',
        workdir: str = "."
) -> None:
    """
    Call gffread to produce cds/protein sequence file (if you have gffread), by running\n
    `abs_gffread_path gff3_file_path -g fna_file_path -x cds_out`\n
    and \n
    `abs_gffread_path gff3_file_path -g fna_file_path -y prot_out`.\n
    After calling, it will remove temporary fasta index file (.fna.fai).
    """
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, '-x', cds_out],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   cwd=workdir)
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, "-y", pep_out],
                   cwd=workdir)
    os.remove(f'{workdir}/{fna_file_path}.fai')


from typing import Literal

class Diamond:
    """
    Diamond is a high-perference program for big batch BLASTp/BLASTx tasks.
    """
    def __init__(
        self,
        fasta_query_path: str,
        fasta_target_path: str,
        diamond_path: str = "diamond",
        verbose: bool = False,
        threads: int = 6,
        min_evalue: float = 1e-3,
        max_target_seqs: int = 25,
        sensitivity: Literal[
            "faster", "fast", "mid-sensitive", "sensitive",
            "more-sensitive", "very-sensitive", "ultra-sensitive"] = "ultra-sensitive",
        outfmt: int = 6,
        no_self_hits: bool = True,
    ) -> None:
        self.tmpdir = ""
        self.makedb_command = [
            diamond_path,
            "makedb",
            f"--threads {threads}",
            f"--db {self.tmpdir}/diamond.db",
            f"--verbose {verbose}",
        ]
        self.diamond_command = [
            diamond_path,
            "blastp",
            f"--threads {threads}",
        ]
