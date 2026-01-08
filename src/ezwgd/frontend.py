'''
Program frontend.
'''

import subprocess
import os
import shutil

from typing import Literal
from pydantic import BaseModel, Field, ConfigDict, field_validator

from ._custom import ProgramNotFoundError


def gffread_cds_pep_extract(
        gff3_file_path: str,
        fna_file_path: str,
        cds_out: str,
        pep_out: str,
        abs_gffread_path: str = 'gffread',
        workdir: str = '.'
) -> None:
    '''
    Call gffread to produce cds/protein sequence file (if you have gffread), by running\n
    `abs_gffread_path gff3_file_path -g fna_file_path -x cds_out`\n
    and \n
    `abs_gffread_path gff3_file_path -g fna_file_path -y prot_out`.\n
    After calling, it will remove temporary fasta index file (.fna.fai).
    '''
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, '-x', cds_out],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   cwd=workdir)
    subprocess.run([abs_gffread_path, gff3_file_path, '-g', fna_file_path, '-y', pep_out],
                   cwd=workdir)
    os.remove(f'{workdir}/{fna_file_path}.fai')


class Diamond:
    '''
    Diamond is a high-perference program for big batch BLASTp/BLASTx tasks.
    '''
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


# ---Wrapper API for PAML package, rewriten by Pydantic, from Biopython.------------------------------
# ---Map for NCBI GenBank Codon Table to PAML ICODE.--------------------------------------------------
NCBI_TO_PAML_ICODE = {
    '1': 0, 'Standard': 0,
    '2': 1, 'Vertebrate Mt': 1,
    '3': 2, 'Yeast Mt': 2,
    '4': 3, 'Mold Mt': 3,
    '5': 4, 'Invertebrate Mt': 4,
    '6': 5, 'Ciliate Nucl': 5,
    '9': 6, 'Echinoderm Mt': 6,
    '10': 7, 'Euplotid Nucl': 7,
    '12': 8, 'Alternative Yeast Nucl': 8,
    '13': 9, 'Ascidian Mt': 9,
    '15': 10, 'Blepharisma Nucl': 10,
}
class CODEML(BaseModel):
    '''
    '''
    # forbid unknown parameters.
    model_config = ConfigDict(extra='forbid')


    noisy: Literal[0, 1, 2, 3 ,9] = Field(default=0,)
    'Level of detail in iteration progress display.'

    verbose: bool = Field(default=False)
    'Whether to display terminal output in detail.'

    runmode: Literal['usertree', 'semiauto', 'fullauto', 'stepwise', 'simpNNI', 'treeNNI',
                     'pairwise', 'pairwiseBayesian'] = Field(default='usertree',)
    'Program operation mode. Note that the program author has stated that "options with NNI may not function as expected."'
    
    t_gamma: tuple[float, float] = Field(default=(1.1, 1.1),)

    omega_gamma: tuple[float, float] = Field(default=(1.1, 1.2))

    seqtype: Literal['codon', 'aa', 'codon2aa','1', '2', '3', 1, 2 ,3] = Field(default='codon')
    '''
    How to use the sequence data.\n
    - `codon`(`1`): these sequences are codons(DNAs).\n
    - `aa`(`2`): these sequences are amino acids.\n
    - `codon2aa`(`3`): these sequences are codons, but will be translated to amino acids.
    '''

    @field_validator('seqtype', mode='before')
    def parse_seqtype(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'codon': '1', 'aa': '2', 'codon2aa': '3'}.get(v)
    
    CodonFreq: Literal['Average', 'F1X4', 'F3X4', 'Fcodon', 'F1X4MG', 'F3X4MG', 'FMutSel0', 'FMutSel',
                       '0', '1', '2', '3', '4', '5', '6', '7', 0, 1, 2, 3, 4, 5, 6, 7] = Field(default='F3X4',)
    'Codon substitution matrix.'

    @field_validator('CodonFreq', mode='before')
    def parse_CodonFreq(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'Average': '0', 'F1X4': '1', 'F3X4': '2', 'Fcodon': '3', 'F1X4MG': '4', 'F3X4MG': '5', 'FMutSel0': '6', 'FMutSel': '7'}.get(v)

    ndata: int = Field(default=1,)
    'Number of datasets.'

    clock: Literal['no', 'global', 'local', '0', '1', '2', 0, 1, 2] = Field(default='no',)
    '''
    Use (or not use) the molecular clock hypothesis, and if you use, which molecular clock has been choosed.
    - `no`(`0`)
    - `global`(`1`)
    - `local`(`2`)
    '''
 
    @field_validator('clock', mode='before')
    def parse_clock(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'no': '0', 'global': '1', 'local': '2'}.get(v)

    model: Literal[
        'one-ratio', 'free-ratios', 'more-per-branch', 'poisson', 'proportional', 'Empirical', 'Empirical+F', 'FromCodon0', 'FromCodon1', 'REVaa0', 'REVaa',
        '0', '1', '2', '3', '5', '6', '8', '9', 0, 1, 2, 3, 5, 6, 8, 9] = Field()
    '''
    Model type.\n
    Models for codons:\n
    0:one, 1:b, 2:2 or more dN/dS ratios for branches
    Models for AAs or codon-translated AAs\n
    0:poisson, 1:proportional,2:Empirical,3:Empirical+F<br/>* 5:FromCodon0, 6:FromCodon1, 8:REVaa_0, 9:REVaa(nr=189)
    '''

    @field_validator('model', mode='before')
    def parse_model(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'one-ratio': '0', 'free-ratios': '1', 'more-per-branch': '2',
                    'poisson': '0', 'proportional': '1', 'Empirical': '2', 'Empirical+F': '3', 'FromCodon0': '5', 'FromCodon1': '6', 'REVaa0': '8', 'REVaa': '9',}.get(v)

    NSsites: Literal['one-omega', 'netural', 'selection', 'discrete', 'freqs', 'gamma', '2gamma', 'beta', 'beta-omega', 'beta-gamma', 'beta-gamma1', 'beta&normal>1', '0&2normal>1', '3normal>0',
                     '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] = Field()
    '''

    '''

    @field_validator('NSsites', mode='before')
    def parse_NSsites(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'one-omega': '0', 'netural': '1', 'selection': '2', 'discrete': '3', 'freqs': '4', 'gamma': '5', '2gamma': '6', 'beta': '7', 'beta-omega': '8', 'beta-gamma': '9',
                    'beta-gamma1': '10', 'beta&normal>1': '11', '0&2normal>1': '12', '3normal>0': '13',}.get(v)

    estFreq: Literal['observed', 'estimate', '0', '1', 0, 1] = Field(default='observed',)
    '''
    Which codon frequency should be used to construct the matrix.\n
    - `observed`(`0`): Using observed frequency.\n
    - `estimate`(`1`): Using ML-estimated frequency.
    '''

    @field_validator('estFreq', mode='before')
    def parse_estFreq(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'observed': '0', 'estimate': '1'}.get(v)

    icode: Literal['Standard', 'Vertebrate Mt', 'Yeast Mt', 'Mold Mt', 'Invertebrate Mt', 'Ciliate Nucl', 'Echinoderm Mt', 'Euplotid Nucl', 'Alternative Yeast Nucl', 'Ascidian Mt', 'Blepharisma Nucl',
                   '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] = Field(default='Standard',)
    'Codon type. You do not need to modify it, unless you have specific requirements.'

    @field_validator('icode', mode='before')
    def parse_icode(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'Standard': '0', 'Vertebrate Mt': '1', 'Yeast Mt': '2', 'Mold Mt': '3', 'Invertebrate Mt': '4', 'Ciliate Nucl': '5',
                    'Echinoderm Mt': '6', 'Euplotid Nucl': '7', 'Alternative Yeast Nucl': '8', 'Ascidian Mt': '9', 'Blepharisma Nucl': '10'}.get(v)

    Mgege: Literal['rates', 'separate'] = 'rates'
    fix_kappa: bool = Field(
        default=False,
        description='',
    )
    kappa: float = Field(
        default=2,
        description='',
    )
    fix_omega: bool = Field(
        default=False,
        description='',
    )
    omega: float = Field(
        default=0.4,
        description='',
    )
    fix_alpha: bool = Field(
        default=True,
        description='',
    )
    alpha: float = Field(default=0,)
    malpha: bool = Field(
        default=True,
        description='',
    )
    fix_rho: bool = Field(
        default=True,
        description='',
    )
    rho: float = Field(
        default=0,
        description='',
    )
    ncatG: int = Field(
        default=3,
        description='',
    )
    getSE: bool = Field(
        default=False,
    )
    RateAncestor: Literal[0, 1, 2] = Field(
        default=0,
        description='',
    )
    Small_Diff: float = Field(
        default=0.5e-6,
        description='',
    )
    cleandata: bool = Field(
        default=False,
        description='',
    )
    fix_blength: Literal['ignore', 'random', 'initial', 'fixed', 'proportional', '0', '-1', '1', '2', '3', 0, -1, 1, 2, 3] = Field(default='ignore',)
    '''
    
    '''

    @field_validator('fix_blength', mode='before')
    def parse_fix_blength(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'ignore': '0', 'random': '-1', 'initial': '1', 'fixed': '2', 'proportional': '3'}.get(v)

    method: Literal['simultaneous', 'eachbranch', '0', '1', 0, 1] = Field(default='simultaneous',)
    '''

    - `simultaneous`(`0`)
    - `eachbranch`(`1`)
    '''

    @field_validator('method', mode='before')
    def parse_method(cls, v):
        if isinstance(v, int):
            return str(v)
        else:
            return {'simultaneous': 0, 'eachbranch': 1}.get(v)


# ----------------------------------------------------------------------------------------------------
class ConfigYN00:
    '''
    Define the parameters used by the yn00 program. If you want custom those parameters, you can instantiate this,
    change you want, and passed as an argument for `evo.dnds.calc_yn00`.\n
    These parameters correspond exactly to those required by the yn00 program, but include additional details.\n
    Furthermore, when running multiple processes, these parameters will be applied to all processes.\n

    The descriptions of each attribute/parameter are as follows:\n
    - `verbose`: (Original) Display detailed output on the terminal. Use this if you want debug.\n
    - `codon_table`: (Original, changed from paml `icode`) Type of codon table, by NCBI GenBank Codon Type index.\n
    - `weighting`: (Original) Weighting based on protein alignment quality.\n
    - `commonf3x4`: (Original) Using the common F3x4 codon frequency.\n
    - `align_method`: Protein alignment methods.\n
    - `detail_result`: More detailed result with kappa, omega and standard error (SE).\n
    - `process_number`: Number of concurrent processes.\n
    - `isaligned`: This pair of protein sequences has been aligned previously.\n
    '''
    def __init__(
        self,
        verbose: bool = False,
        codon_table: str = '1',
        weighting: bool = False,
        commonf3x4: bool = False,
        align_method: Literal['muscle', 'mafft'] = 'muscle',
        detail_result: bool = False,
        processes_number: int = 6,
        isaligned: bool = False,
        ) -> None:
        self.verbose = int(verbose),
        self.icode = NCBI_TO_PAML_ICODE.get(str(codon_table), 0),
        self.weighting = int(weighting),
        self.commonf3x4 = int(commonf3x4)
        self.align_method: Literal['muscle', 'mafft'] = align_method
        self.detail_result = detail_result
        self.processes_number = processes_number
        self.isaligned = isaligned

        # Fast check (for executable files).
        for program in ['yn00', align_method]:
            if shutil.which(program) is None:
                raise ProgramNotFoundError(program)


class Chi2:
    pass
    # TODO chi2