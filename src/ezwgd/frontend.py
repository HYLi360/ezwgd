#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import subprocess
import os
import shutil

from dataclasses import dataclass, fields
from tempfile import TemporaryDirectory

from typing import Literal

from ._custom import ProgramNotFoundError

from .utils.parse import codeml_pairwise


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


# ---Wrapper API for PAML package, rewriten by Pydantic, from Biopython.------------------------------
# ---Map for NCBI GenBank Codon Table to PAML ICODE.--------------------------------------------------
NCBI_TO_PAML_ICODE = {
    '1': '0', 'Standard': '0',
    '2': '1', 'Vertebrate Mt': '1',
    '3': '2', 'Yeast Mt': '2',
    '4': '3', 'Mold Mt': '3',
    '5': '4', 'Invertebrate Mt': '4',
    '6': '5', 'Ciliate Nucl': '5',
    '9': '6', 'Echinoderm Mt': '6',
    '10': '7', 'Euplotid Nucl': '7',
    '12': '8', 'Alternative Yeast Nucl': '8',
    '13': '9', 'Ascidian Mt': '9',
    '15': '10', 'Blepharisma Nucl': '10',
}

_NOISY_LITERAL = Literal['none', 'low', 'medium', 'high', 'max']
_NOISY_MAP = {'none': 0, 'low': 1, 'medium': 2, 'high': 3, 'max': 9}
_RUNMODE_LITERAL = Literal['usertree', 'semiauto', 'fullauto', 'stepwise', 'simpNNI', 'treeNNI', 'pairwise']
_RUNMODE_MAP = {'usertree': 0, 'semiauto': 1, 'fullauto': 2, 'stepwise': 3, 'simpNNI': 4, 'treeNNI': 5, 'pairwise': -2}
_SEQTYPE_LITERAL = Literal['codon', 'aa', 'codon2aa']
_SEQTYPE_MAP = {'codon': 1, 'aa': 2, 'codon2aa': 3}
_CODONFREQ_LITERAL = Literal['Average', 'F1X4', 'F3X4', 'Fcodon', 'F1X4MG', 'F3X4MG', 'FMutSel0', 'FMutSel']
_CODONFREQ_MAP = {'Average': 0, 'F1X4': 1, 'F3X4': 2, 'Fcodon': 3, 'F1X4MG': 4, 'F3X4MG': 5, 'FMutSel0': 6,
                  'FMutSel': 7}
_CLOCK_LITERAL = Literal['no', 'global', 'local']
_CLOCK_MAP = {'no': 0, 'global': 1, 'local': 2}
_MODEL_LITERAL = Literal[
    'one-ratio', 'free-ratios', 'more-per-branch', 'poisson', 'proportional', 'Empirical', 'Empirical+F', 'FromCodon0', 'FromCodon1', 'REVaa0', 'REVaa']
_MODEL_MAP = {'one-ratio': 0, 'free-ratios': 1, 'more-per-branch': 2,
              'poisson': 0, 'proportional': 1, 'Empirical': 2, 'Empirical+F': 3, 'FromCodon0': 5, 'FromCodon1': 6,
              'REVaa0': 8, 'REVaa': 9}
_NSSITES_LITERAL = Literal[
    'one-omega', 'netural', 'selection', 'discrete', 'freqs', 'gamma', '2gamma', 'beta', 'beta-omega', 'beta-gamma', 'beta-gamma1', 'beta-normal>1', '0-2normal>1', '3normal>0']
_NSSITES_MAP = {'one-omega': 0, 'netural': 1, 'selection': 2, 'discrete': 3, 'freqs': 4, 'gamma': 5, '2gamma': 6,
                'beta': 7, 'beta-omega': 8, 'beta-gamma': 9, 'beta-gamma1': 10, 'beta-normal>1': 11, '0-2normal>1': 12,
                '3normal>0': 13}
_ESTFREQ_LITERAL = Literal['observed', 'estimate']
_ESTFREQ_MAP = {'observed': 0, 'estimate': 1}
_ICODE_LITERAL = Literal[
    'Standard', 'Vertebrate Mt', 'Yeast Mt', 'Mold Mt', 'Invertebrate Mt', 'Ciliate Nucl', 'Echinoderm Mt', 'Euplotid Nucl', 'Alternative Yeast Nucl', 'Ascidian Mt', 'Blepharisma Nucl']
_ICODE_MAP = {'Standard': 0, 'Vertebrate Mt': 1, 'Yeast Mt': 2, 'Mold Mt': 3, 'Invertebrate Mt': 4, 'Ciliate Nucl': 5,
              'Echinoderm Mt': 6, 'Euplotid Nucl': 7, 'Alternative Yeast Nucl': 8, 'Ascidian Mt': 9,
              'Blepharisma Nucl': 10}
_MGENE_LITERAL = Literal['rates', 'separate']
_MGENE_MAP = {'rates': 0, 'separate': 1}
_RATEANCESTOR_LITERAL = Literal[0, 1, 2]
_RATEANCESTOR_MAP = {0: 0, 1: 1, 2: 2}
_FIX_BLENGTH_LITERAL = Literal['ignore', 'random', 'initial', 'fixed', 'proportional']
_FIX_BLENGTH_MAP = {'ignore': 0, 'random': -1, 'initial': 1, 'fixed': 2, 'proportional': 3}
_METHOD_LITERAL = Literal['simultaneous', 'eachbranch']
_METHOD_MAP = {'simultaneous': 0, 'eachbranch': 1}


@dataclass(kw_only=True, )
class ConfigCODEML:
    seqfile: str = ''
    outfile: str = ''

    noisy: _NOISY_LITERAL = 'none'
    '''Level of detail in iteration progress display.'''

    verbose: bool = True
    '''Whether to display terminal output in detail.'''

    runmode: _RUNMODE_LITERAL = 'usertree'
    '''Program operation mode. Note that the program author has stated that "options with NNI may not function as expected."'''

    seqtype: _SEQTYPE_LITERAL = 'codon'
    '''
    How to use the sequence data.\n
    - `codon`: these sequences are codons(DNAs).\n
    - `aa`: these sequences are amino acids.\n
    - `codon2aa`: these sequences are codons, but will be translated to amino acids.
    '''

    CodonFreq: _CODONFREQ_LITERAL = 'F3X4'
    '''Codon substitution matrix.'''

    ndata: int = 1
    '''Number of datasets.'''

    clock: _CLOCK_LITERAL = 'no'
    '''Use (or not use) the molecular clock hypothesis, and if you use, which molecular clock has been choosed.'''

    model: _MODEL_LITERAL = 'one-omega'
    '''
    Model type.\n
    Models for codons:\n`one-ratio`, `free-ratios`, `more-per-branch`\n
    Models for AAs or codon-translated AAs:\n`poisson`, `proportional`, `Empirical`, `Empirical+F`, `FromCodon0`, `FromCodon1`, `REVaa0`, `REVaa`
    '''

    NSsites: _NSSITES_LITERAL = 'one-omega'
    ''''''

    estFreq: _ESTFREQ_LITERAL = 'observed'
    '''Which codon frequency should be used to construct the matrix.'''

    icode: _ICODE_LITERAL = 'Standard'
    '''Codon type. You do not need to modify it, unless you have specific requirements.'''

    Mgene: _MGENE_LITERAL = 'rates'
    ''''''

    fix_kappa: bool = False
    ''''''

    kappa: float = 2
    ''''''

    fix_omega: bool = False
    ''''''

    omega: float = 0.4
    ''''''

    fix_alpha: bool = True
    ''''''

    alpha: float = 0
    ''''''

    Malpha: bool = True
    ''''''

    fix_rho: bool = True
    ''''''

    rho: float = 0
    ''''''

    ncatG: int = 3
    ''''''

    getSE: bool = False
    ''''''

    RateAncestor: _RATEANCESTOR_LITERAL = 0
    ''''''

    Small_Diff: float = 5e-7
    ''''''

    cleandata: bool = True
    ''''''

    fix_blength: _FIX_BLENGTH_LITERAL = 'ignore'
    ''''''

    method: _METHOD_LITERAL = 'simultaneous'
    ''''''


@dataclass(kw_only=True, )
class CODEML:
    """
    This class translates the `ConfigCODEML` property of semantic versioning into a string,
    facilitating its output as a control file.\n
    Advanced users, who accustomed to using `CODEML` in traditional ways, can also directly use this class.
    """
    seqfile: str = ''
    outfile: str = ''
    noisy: int = 0
    verbose: int = 0
    runmode: int = 0
    seqtype: int = 0
    CodonFreq: int = 2
    ndata: int = 1
    clock: int = 0
    model: int = 0
    NSsites: int = 0
    estFreq: int = 0
    icode: int = 0
    Mgene: int = 0
    fix_kappa: int = 0
    kappa: float = 2
    fix_omega: int = 0
    omega: float = 0.4
    fix_alpha: int = 1
    alpha: float = 0
    Malpha: int = 1
    fix_rho: int = 1
    rho: float = 0
    ncatG: int = 3
    getSE: float = 0
    RateAncestor: float = 0
    Small_Diff: float = 5e-7
    cleandata: int = 1
    fix_blength: int = 0
    method: int = 0

    @staticmethod
    def _str2int(param_value, mapper: dict,) -> int:
        return mapper.get(param_value, 999)

    @staticmethod
    def _bool2int(param_value: bool,) -> int:
        return int(param_value)

    @classmethod
    def escape(cls, config: ConfigCODEML) -> "CODEML":
        """
        Use this if you want to escape (from `ConfigCODEML` to `CODEML`):\n
        `adv_cfg = CODEML.escape(cfg)`\n
        which the `cfg` is a instance of `ConfigCODEML`.
        """
        return cls(
            seqfile=config.seqfile,
            outfile=config.outfile,
            noisy=cls._str2int(config.noisy, _NOISY_MAP),
            verbose=cls._bool2int(config.verbose),
            runmode=cls._str2int(config.runmode, _RUNMODE_MAP),
            seqtype=cls._str2int(config.seqtype, _SEQTYPE_MAP),
            CodonFreq=cls._str2int(config.CodonFreq, _CODONFREQ_MAP),
            ndata=config.ndata,
            clock=cls._str2int(config.clock, _CLOCK_MAP),
            model=cls._str2int(config.model, _MODEL_MAP),
            NSsites=cls._str2int(config.NSsites, _NSSITES_MAP),
            estFreq=cls._str2int(config.estFreq, _ESTFREQ_MAP),
            icode=cls._str2int(config.icode, _ICODE_MAP),
            Mgene=cls._str2int(config.Mgene, _MGENE_MAP),
            fix_kappa=cls._bool2int(config.fix_kappa),
            kappa=config.kappa,
            fix_omega=cls._bool2int(config.fix_omega),
            omega=config.omega,
            fix_alpha=cls._bool2int(config.fix_alpha),
            alpha=config.alpha,
            Malpha=cls._bool2int(config.Malpha),
            fix_rho=cls._bool2int(config.fix_rho),
            rho=config.rho,
            ncatG=config.ncatG,
            getSE=cls._bool2int(config.getSE),
            RateAncestor=cls._str2int(config.RateAncestor, _RATEANCESTOR_MAP),
            Small_Diff=config.Small_Diff,
            cleandata=cls._bool2int(config.cleandata),
            fix_blength=cls._str2int(config.fix_blength, _FIX_BLENGTH_MAP),
            method=cls._str2int(config.method, _METHOD_MAP),
        )

    def __str__(self) -> str:
        lines = []
        for f in fields(self):
            lines.append(f"{f.name} = {getattr(self, f.name)}")
        return "\n".join(lines) + "\n"

    def preset_dnds(self,
                    work_dir: str,
                    codon_data: str,
                    codeml_cmd: str = 'codeml'
                    ) -> dict:
        """
        A preset to calculate dN/dS.\n
        Preset parameters:\n
        - `runmode`: pairwise(-2), which DOESN'T need tree file.\n
        - `seqtype`: codon(1).\n
        - `model`: free-ratios(1).\n
        Those parameters will apply directly and cannot change, but doesn't affect others.\n
        `{'N': xxx, 'S': xxx, 'dN': xxx, 'dS': xxx, 'omega': xxx, 't': xxx}`
        Args:
            codon_data (str): codon sequence data (after codon align).
            work_dir (str): work directory.
            codeml_cmd (str): The command (or path) of codeml.
        """
        # preset
        self.seqfile = 'codon.aln'
        self.outfile = 'res.txt'
        self.runmode = -2
        self.seqtype = 1
        self.model = 1
        # write codon sequence file
        with open(f'{work_dir}/{self.seqfile}', 'w') as file: file.write(codon_data)
        # write control
        with open(f'{work_dir}/codeml.ctl', 'w') as file: file.write(str(self))
        # run codeml
        subprocess.run(
            args=[codeml_cmd],
            cwd=work_dir,
            stdout=subprocess.DEVNULL if not self.verbose else None
        )
        # read result
        return codeml_pairwise(f'{work_dir}/rst')


# ----------------------------------------------------------------------------------------------------
class ConfigYN00:
    """
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
    """

    def __init__(
            self,
            verbose: bool = False,
            codon_table: str = '1',
            weighting: bool = False,
            commonf3x4: bool = False,
            align_method: Literal['muscle', 'mafft', None] = 'muscle',
            processes_number: int = 6,
    ) -> None:
        """

        Args:
            verbose ():
            codon_table ():
            weighting ():
            commonf3x4 ():
            align_method ():
            processes_number ():
        """
        self.verbose = int(verbose),
        self.icode = NCBI_TO_PAML_ICODE.get(str(codon_table), 0),
        self.weighting = int(weighting),
        self.commonf3x4 = int(commonf3x4)
        self.align_method: Literal['muscle', 'mafft', None] = align_method
        self.processes_number = processes_number

        # Fast check (for executable files).
        for program in ['yn00', align_method]:
            if shutil.which(program) is None:
                raise ProgramNotFoundError(program)


class Chi2:
    pass
    # TODO chi2
