# ---Import.------------------------------------------------------------------------------------------
from multiprocessing import Pool
import tempfile

from typing import Literal
from math import ceil

from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from Bio.SeqRecord import SeqRecord

from ..frontend import muscle5, mafft, ConfigYN00


# ----------------
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


# ---Dispatch by algorithm/executable program name (MUSCLE/MAFFT).------------------------------------
def _alignp(
        method: Literal["muscle", "mafft"],
        workdir: str,
        pep_file_name: str,
        ) -> tuple[SeqRecord, SeqRecord]:
    match method:
        case "muscle":
            muscle5(workdir, pep_file_name, "align.faa")
        case "mafft":
            mafft(workdir, pep_file_name, "align.faa")
    return tuple(SeqIO.parse(f"{workdir}/align.faa", "fasta"))


# ---calculate by ng86 and yn00.----------------------------------------------------------------------
def _run_yn00(
        config: ConfigYN00,
        workdir: str
        ) -> dict:
    proc_yn00 = yn00.Yn00(
        working_dir = workdir,
        alignment = f"{workdir}/codon.aln",
        out_file = f"{workdir}/res.txt",)

    proc_yn00.set_options(
        verbose = config.verbose,
        icode = config.icode,
        weighting = config.weighting,
        commonf3x4 = config.commonf3x4)

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
            prot_aligned = _alignp(config.align_method, tmpdir, "prot.faa")

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
        config: ConfigYN00 = ConfigYN00(),
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

