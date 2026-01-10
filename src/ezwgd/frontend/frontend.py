#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import subprocess
import os


from typing import Literal




def gffread_cds_pep_extract(
        gff3_file_path: str,
        fna_file_path: str,
        cds_out: str,
        pep_out: str,
        abs_gffread_path: str = 'gffread',
        workdir: str = '.'
) -> None:
    """
    Call `gffread` to produce cds/protein sequence file (if you have gffread), by running\n
    `abs_gffread_path gff3_file_path -g fna_file_path -x cds_out`\n
    and \n
    `abs_gffread_path gff3_file_path -g fna_file_path -y prot_out`.\n
    After calling, it will remove temporary fasta index file (.fna.fai).
    Args:
        gff3_file_path (`str`):
        fna_file_path (`str`):
        cds_out (`str`):
        pep_out (`str`):
        abs_gffread_path (`str`):
        workdir (`str`):
    """
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, '-x', cds_out],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   cwd=workdir)
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, '-y', pep_out],
                   cwd=workdir)
    os.remove(f'{workdir}/{fna_file_path}.fai')


class Diamond:
    def __init__(
            self,
            fasta_query_path: str,
            fasta_target_path: str,
            diamond_path: str = 'diamond',
            verbose: bool = False,
            threads: int = 6,
            min_evalue: float = 1e-3,
            max_target_seqs: int = 25,
            sensitivity: Literal[
                'faster', 'fast', 'mid-sensitive', 'sensitive',
                'more-sensitive', 'very-sensitive', 'ultra-sensitive'] = 'ultra-sensitive',
            outfmt: int = 6,
            no_self_hits: bool = True,
    ) -> None:
        self.tmpdir = ''
        self.makedb_command = [
            diamond_path,
            'makedb',
            f'--threads {threads}',
            f'--db {self.tmpdir}/diamond.db',
            f'--verbose {verbose}',
        ]
        self.diamond_command = [
            diamond_path,
            'blastp',
            f'--threads {threads}',
        ]


def muscle5(
        workdir: str,
        pep_file_name: str,
        res_file_name: str
) -> None:
    subprocess.run(
        ['muscle', '-align', pep_file_name, '-output', res_file_name],
        stderr=subprocess.PIPE,
        text=True,
        cwd=workdir,
    )


def mafft(
        workdir: str,
        pep_file_name: str,
        res_file_name: str
) -> None:
    with open(f'{workdir}/{res_file_name}', 'w') as outfile:
        subprocess.run(
            ['mafft', '--auto', pep_file_name],
            stdout=outfile,
            stderr=subprocess.PIPE,
            text=True,
            cwd=workdir,
        )

