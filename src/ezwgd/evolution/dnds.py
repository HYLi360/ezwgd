# ---Import.------------------------------------------------------------------------------------------
from multiprocessing import Pool
import subprocess
import tempfile
import shutil

from typing import Literal
from math import ceil

from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from Bio.SeqRecord import SeqRecord

from .._custom import ProgramNotFoundError

# ---Map for NCBI GenBank Codon Table to PAML ICODE.--------------------------------------------------
NCBI_TO_PAML_ICODE = {
    1: 0,
    "Standard": 0,
    2: 1,
    "Vertebrate Mt": 1,
    3: 2,
    "Yeast Mt": 2,
    4: 3,
    "Mold Mt": 3,
    5: 4,
    "Invertebrate Mt": 4,
    6: 5,
    "Ciliate Nucl": 5,
    9: 6,
    "Echinoderm Mt": 6,
    10: 7,
    "Euplotid Nucl": 7,
    12: 8,
    "Alternative Yeast Nucl": 8,
    13: 9,
    "Ascidian Mt": 9,
    15: 10,
    "Blepharisma Nucl": 10,
}


# ----------------------------------------------------------------------------------------------------
class ConfigYN00:
    """
    Define the parameters used by the yn00 program. Before executing cacl_yn00, the class must be instantiated and passed as an argument.\n
    These parameters correspond exactly to those required by the yn00 program, but include additional details.\n
    Furthermore, when running multiple processes, these parameters will be applied to all processes.\n
    If no adjustments are made, it will use default parameters.\n

    The descriptions of each attribute/parameter are as follows:\n
    - `verbose`: (Original) Display detailed output on the terminal.\n
    - `icode`: (Original) Type of codon table.\n
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
        codon_table: str = "1",
        weighting: bool = False,
        commonf3x4: bool = False,
        align_method: Literal["muscle", "mafft"] = "muscle",
        detail_result: bool = False,
        processes_number: int = 6,
        isaligned: bool = False,
    ) -> None:
        self.verbose = verbose
        self.codon_table = codon_table
        self.paml_icode = NCBI_TO_PAML_ICODE.get(str(codon_table), 0)
        self.weighting = weighting
        self.commonf3x4 = commonf3x4
        self.align_method: Literal["muscle", "mafft"] = align_method
        self.detail_result = detail_result
        self.processes_number = processes_number
        self.isaligned = isaligned

        # Fast check (for executable files).
        for program in ["yn00", align_method]:
            if shutil.which(program) is None:
                raise ProgramNotFoundError(program)

def _codon_align(
        cds_records: tuple[SeqRecord, SeqRecord],
        protein_aligned: tuple[SeqRecord, SeqRecord],
        workdir: str,
        ) -> None:
        # init
        cds1, cds2 = str(cds_records[0].seq), str(cds_records[1].seq)
        prot1, prot2 = str(protein_aligned[0].seq), str(protein_aligned[1].seq),
        cdsp1, cdsp2 = 0, 0
        codons_buf1, codons_buf2 = [], []

        # len(prot1) == len(prot2)
        for i in range(len(prot1)):
            prot_site1, prot_site2 = prot1[i], prot2[i]
            # Extract two codons when...
            # prot1 ......A......
            # prot2 ......A......
            if (prot_site1 != "-") and (prot_site2 != "-"):
                # pick codons
                codon1, codon2 = (
                    cds1[cdsp1 : cdsp1 + 3],
                    cds2[cdsp2 : cdsp2 + 3],
                )
                # concat codons
                codons_buf1.append(codon1)
                codons_buf2.append(codon2)
                # move the pointers
                cdsp1 += 3
                cdsp2 += 3
            # Skip when...
            # prot1 ......-......
            # prot2 ......A......
            elif prot_site1 != "-":
                cdsp1 += 3
            elif prot_site2 != "-":
                cdsp2 += 3
            continue

        # generate paml input
        codons_buf1 = "".join(codons_buf1)
        codons_buf2 = "".join(codons_buf2)
        with open(f"{workdir}/codon.aln", "w+") as f:
            # Structure of paml input: header (sequence count (4-char length) + paired base count (7-char length)
            # seq1
            # (Detailed sequence)
            # seq2
            # (Detailed sequence)
            # Concat with "\n" and join method.
            header = f"{2: 4}{len(codons_buf1): 7}"
            f.write("\n".join([header, "seq1", codons_buf1, "seq2", codons_buf2]))


# ---PIs for various protein alignment methods.-------------------------------------------------------
def _muscle5(workdir: str):
    subprocess.run(
        ["muscle", "-align", "prot.faa", "-output", "align.faa"],
        stderr=subprocess.PIPE,
        text=True,
        cwd=workdir,
    )


def _mafft(workdir: str):
    with open(f"{workdir}/align.faa", "w") as outfile:
        subprocess.run(
            ["mafft", "--auto", "prot.faa"],
            stdout=outfile,
            stderr=subprocess.PIPE,
            text=True,
            cwd=workdir,
        )


# ---Dispatch by algorithm/executable program name (MUSCLE/MAFFT).------------------------------------
def _alignp(
    method: Literal["muscle", "mafft"], workdir: str
) -> tuple[SeqRecord, SeqRecord]:
    match method:
        case "muscle":
            _muscle5(workdir)
        case "mafft":
            _mafft(workdir)
    return tuple(SeqIO.parse(f"{workdir}/align.faa", "fasta"))


# ---calculate ng86 and yn00.-------------------------------------------------------------------------
def _run_yn00(config: ConfigYN00, workdir: str) -> dict:
    proc_yn00 = yn00.Yn00(
        working_dir = workdir,
        alignment = f"{workdir}/codon.aln",
        out_file = f"{workdir}/res.txt",
    )

    proc_yn00.set_options(
        verbose = config.verbose,
        icode = config.paml_icode,
        weighting = config.weighting,
        commonf3x4 = config.commonf3x4,
    )

    res_total = proc_yn00.run(command="yn00", parse=True)
    # Block result if returns None
    if res_total is None:
        raise RuntimeError(
            "yn00 failed to produce result. Check input sequences for internal stop codons."
        )

    res = {
        "NG86": res_total["seq1"]["seq2"]["NG86"],
        "YN00": res_total["seq1"]["seq2"]["YN00"],
    }
    return (
        res
        if config.detail_result
        else {
            "NG86": {"dN": res["NG86"]["dN"], "dS": res["NG86"]["dS"]},
            "YN00": {"dN": res["YN00"]["dN"], "dS": res["YN00"]["dS"]},
        }
    )


# ----------------------------------------------------------------------------------------------------
def _calc_yn00_tasks(
    cds: tuple[SeqRecord, SeqRecord],
    prot: tuple[SeqRecord, SeqRecord],
    config: ConfigYN00,
) -> dict:
    # Tasks should focus on executing the task at hand, rather than preparing in advance to handle exceptions.
    # Create a temporary folder for single tasks.
    with tempfile.TemporaryDirectory(dir="/dev/shm") as tmpdir:
        # Step0: Write records as file.
        SeqIO.write(cds, f"{tmpdir}/cds.fna", "fasta")
        if config.isaligned:
            prot_aligned = prot
        else:
            # Step1: Protein Alignment (if needed).
            SeqIO.write(prot, f"{tmpdir}/prot.faa", "fasta")
            prot_aligned = _alignp(config.align_method, tmpdir)

        # Step2: Codon Alignment.
        _codon_align(cds, prot_aligned, tmpdir)

        # Step3: do yn00.
        result = _run_yn00(config, tmpdir)
        result["seqid"] = [cds[0].id, cds[1].id]

        return result

def calc_yn00(
    cds1: str | list[SeqRecord],
    cds2: str | list[SeqRecord],
    prot1: str | list[SeqRecord],
    prot2: str | list[SeqRecord],
    config: ConfigYN00,
    ):
    """
    Calculate the dN and dS values between homologous gene pairs by yn00.

    Args:
        cds:
        prot:
        config:
    
    Return:
        
    """
    total_records = []
    for records in [cds1, cds2, prot1, prot2]:
        total_records.append([record for record in SeqIO.parse(records, "fasta")] if isinstance(records, str) else records)

    cds, prot = list(zip(total_records[0], total_records[1])), list(zip(total_records[2], total_records[3]))

    with Pool(processes=config.processes_number) as pool:
        args_iterable = ((cds[i], prot[i], config) for i in range(len(total_records[0])))
        results = pool.starmap(_calc_yn00_tasks, args_iterable, chunksize=ceil(len(cds)/config.processes_number))
        return results
