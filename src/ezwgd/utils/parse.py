from ezwgd._custom import (
    SeqFileNotFoundError,
    GFF3FileNotFoundError,
)

import ezwgd
from ezwgd.utils import pairwise_re, mcscan_info_line, yn00_res_re

import time
import tempfile
import pandas as pd

from typing import Iterable, Optional
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord

from rich.progress import track

class Genome:
    "write something"
    def __init__(
            self,
            fna_file_path: str,
            gff3_file_path: str
            ) -> None:
        ezwgd.console.log("Packup Genome DNA and GFF3 file.")
        init_start = time.time()
        with ezwgd.console.status("Initialing..."):
            try:
                self.genome_seq = SeqIO.index(fna_file_path, "fasta")
            except FileNotFoundError:
                raise(SeqFileNotFoundError(fna_file_path))

            # Read and format original GFF3.
            try:
                gff = pd.read_csv(
                    gff3_file_path, sep="\t", header=None, comment="#",
                    names=["chr", "source", "type", "start", "end", "score", "strand", "phase", "attribute"],
                    dtype={"chr":"str", "type":"str", "start":"int", "end":"int", "strand":"str", "phase":"str", "attribute":"str"}
                )
            except FileNotFoundError:
                raise(GFF3FileNotFoundError(gff3_file_path))

            gff = gff.sort_values(['chr', 'start'])

            # Filtering useful features, and extract data from original "attribute" columns.
            gff = gff[gff["type"].isin(["gene", "mRNA", "transcript", "CDS"])].copy()
            gff["id"] = gff["attribute"].str.extract(r"ID=([^;]+)")
            gff["parent"] = gff["attribute"].str.extract(r"Parent=([^;]+)")

            # Save for next step.
            self.gff = gff.reset_index(drop=True)

        ezwgd.console.log(f"Initialization finished in {time.time() - init_start:3f} seconds.")

        self._build_indices()

    def _build_indices(self):
        ezwgd.console.log("Rebuild GFF3 tree structure.")
        # Rebuild tree structure (gene->transcript->CDS).
        gff = self.gff
        gff_parse_start = time.time()
        with ezwgd.console.status("Rebuilding..."):
            # -----------------------------------------------
            # CDS table (tx -> CDS)
            cds = gff[gff["type"] == "CDS"].copy()
            cds = cds.rename(columns={"parent":"tx_id"})
            cds = cds[cds["tx_id"].notna()]

            # Parent in CDS may be multiple: Parent=a,b,c
            # explode it
            cds["tx_id"] = cds["tx_id"].str.split(",")
            cds = cds.explode("tx_id", ignore_index=True)
            cds["tx_id"] = cds["tx_id"].astype("str")

            # clean phase
            cds["phase"] = cds["phase"].replace(".", "0").fillna("0").astype(int)

            # store segments per transcript: list of tuples
            cds = cds[["tx_id", "chr", "start", "end", "strand", "phase"]].copy()

            # sorting rule: sort by genomic coordinate start
            cds = cds.sort_values(["chr", "start", "end"]).reset_index(drop=True)

            # -----------------------------------------------
            # transcript(tx) table (gene -> tx)
            tx = gff[gff["type"].isin(["mRNA","transcript"])].copy()
            tx = tx.rename(columns={"id":"tx_id", "parent":"gene_id"})
            tx = tx[tx["tx_id"].notna() & tx["gene_id"].notna()].loc[:,["gene_id", "tx_id"]]

            # -----------------------------------------------
            # total_table. like this:
            #     gene_id    tx_id    chr     start       end strand  phase
            # 0     gene1     rna1      1    371878    371957      +      0
            total = pd.merge(tx, cds, how="inner", on="tx_id").sort_values(["chr", "start"]).reset_index(drop=True)

            # gene -> tx dict
            self._gene_txs = total.groupby("gene_id")["tx_id"].apply(list).to_dict()

            # tx_cds: tx: [(tx, chr, start, end, strand, phase), ......]
            tx_cds_dict = defaultdict(list)
            for row in cds.itertuples(index=False, name=None):
                # row[0] is tx_id
                tx_cds_dict[row[0]].append(row)
            self._tx_cds = dict(tx_cds_dict)

            # -----------------------------------------------
            # gene table (for rename gene and counting genes)
            gene = gff[gff["type"] == "gene"].copy()
            gene = gene.loc[:,["id", "chr", "start", "end", "strand"]].rename(columns={"id": "gene_id"}).sort_values(["chr", "start"]).reset_index(drop=True)
            gene['order'] = gene.groupby('chr').cumcount() + 1

            # coding gene list (as default list to extract)
            genelist = total["gene_id"].drop_duplicates()
            self.genelist = genelist.to_list()

            self.simp_gff = pd.merge(genelist, gene, on="gene_id").loc[:,['gene_id', 'chr', 'order', 'start', 'end', 'strand']]

        ezwgd.console.log(f"Rebuild finished in {time.time() - gff_parse_start:3f} seconds. "
                    f"coding gene entries: {len(self.genelist)}, all gene entries: {len(set(gene['gene_id']))}, "
                    f"transcript entries: {tx.shape[0]}, cds entries: {cds.shape[0]}.")

    # ------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------
    def get_cds(self, tx_id) -> SeqRecord:
        """
        Return spliced CDS sequence for transcript tx_id.
        Handles strand and phase.
        """
        if tx_id not in self._tx_cds:
            raise KeyError(f"Transcript {tx_id} has no CDS features")
        
        #           0    1      2    3       4      5
        # list[(tx_id, chr, start, end, strand, phase), ......]
        segs = self._tx_cds[tx_id]

        # pre-test
        chrom, strand = segs[0][1], segs[0][4]

        if segs[0][1] not in self.genome_seq:
            raise KeyError(f"Chromosome {chrom} not found in fasta")
        
        record_seq = self.genome_seq[chrom]

        if record_seq is None:
            raise KeyError(f"No chromosome sequence for {chrom}")

        # extract start
        pieces = []
        for _, _, start, end, _, _ in segs:
            # GFF is 1-based inclusive; Python is 0-based half-open
            frag = record_seq[start-1:end].seq  
            pieces.append(frag)

        cds_seq = Seq("").join(pieces)

        # handle negative strand
        if strand == "-":
            cds_seq = cds_seq.reverse_complement()

        cds = SeqRecord(seq=cds_seq, id=tx_id, name=tx_id, description="")

        return cds

    def get_pep(self, tx_id, to_stop=True, strict=False) -> SeqRecord:
        """
        Translate CDS to peptide using Biopython.
        strict=True uses cds=True for strict checks.
        """
        cds_seq = self.get_cds(tx_id)
        return self.translate(cds_seq, to_stop, strict)
    
    def translate(self, cds_seq: SeqRecord, to_stop=True, strict=False) -> SeqRecord:
        # If not divisible by 3, you can choose to trim or not.
        # Here is a gentle approach: trim tail to multiple of 3 when not strict.
        if not strict:
            # Defence for not standard CDS length.
            trim = len(cds_seq) % 3
            if trim:
                cds_seq = cds_seq[:-trim]

        try:
            pep = cds_seq.translate(to_stop=to_stop, cds=strict)
            pep.id = cds_seq.id
            pep.name = cds_seq.name
            pep.description = ""
        except Exception as e:
            raise ValueError(f"Translation failed for {cds_seq.id}: {e}")
        return pep

    # ------------------------------------------------------------
    # Longest CDS/Prot per gene
    # ------------------------------------------------------------
    def get_longest(self, gene_id, to_stop=True, strict=False):
        """
        Return (best_tx_id, best_cds_seq, best_pep_seq)
        best determined by CDS nucleotide length (after phase trimming/splicing)
        """
        if gene_id not in self._gene_txs:
            raise KeyError(f"Gene {gene_id} has no transcripts")

        best = None
        for tx_id in self._gene_txs[gene_id]:
            if tx_id not in self._tx_cds:
                continue
            cds = self.get_cds(tx_id)
            if best is None or len(cds) > len(best[1]):
                best = (tx_id, cds)

        if best is None:
            raise ValueError(f"Gene {gene_id} has no CDS")

        tx_id, cds = best
        pep = self.get_pep(tx_id, to_stop=to_stop, strict=strict)
        return tx_id, cds, pep

    def iter_longest_pep(self, to_stop=True, strict=False):
        for gene_id in self._gene_txs:
            try:
                tx_id, cds, pep = self.get_longest(
                    gene_id, to_stop=to_stop, strict=strict
                )
                yield gene_id, tx_id, pep
            except Exception:
                continue
    
    def get_cds_pep(
            self,
            out_cds_path: Optional[str] = None,
            out_pep_path: Optional[str] = None,
            genelist: Optional[Iterable[str]] = None,
            to_stop: bool = True,
            strict: bool = False,
            warp: int = 70,
            ):
        ezwgd.console.log("Extract all CDS and Protein sequences.")
        extract_start = time.time()
        cds_ls, pep_ls = [], []
        genelist = self.genelist if genelist is None else genelist
        for gene_id in track(self.genelist, "Extracting..."):
            tx_id, cds, pep = self.get_longest(gene_id, to_stop, strict)
            cds.id, cds.name, cds.description = gene_id, gene_id, f"tx_id={tx_id}"
            pep.id, pep.name, pep.description = gene_id, gene_id, f"tx_id={tx_id}"
            cds_ls.append(cds)
            pep_ls.append(pep)
        if (out_cds_path is not None) and (out_pep_path is not None):
            cds_handle = FastaWriter(out_cds_path, wrap=warp)
            pep_handle = FastaWriter(out_pep_path, wrap=warp)
            cds_handle.write_file(cds_ls)
            pep_handle.write_file(pep_ls)
        else:
            ezwgd.console.log("You didn't specify the cds & pep output file path. skip writting.")

        ezwgd.console.log(f"Extract finished in {time.time() - extract_start:3f} seconds. "
                    f"{len(cds_ls)} CDS sequences, "
                    f"{len(pep_ls)} Protein sequences exported.")
        
        return cds_ls, pep_ls


class GenomeLight:
    """
    A lighten class to process cds, pep and gff3 file.\n
    It's similar to class `Genome`, but you don't need to use translate tool, such as `gffread`, or `translate` from `Bio`,\n
    so it has much lesser pressure of your hard disk (for I/O and genome sequence storaging usage.).\n

    There have 2 situation for your cds & pep data:\n
    1. These sequences was named by genes' name, not transcripts' name. In this time, we will only make a simple `pd.DataFrame` (`simp_gff`),
    like this:\n
    ```
		gene_id	tx_id	chr	start	end	strand	phase
	0	gene1	rna1	1	371878	371957	+	0
    ```
    2. These sequences was named by transcripts' name. since as one gene can associates with many transcripts, we still have to make a index
    for genes and transcripts, like in `Genome`, and note the gene name as descriptions for every records.\n
    You can change that mode easily by parameter `mode`.
    """
    def __init__(
            self,
            cds_file_path: str,
            pep_file_path: str,
            gff3_file_path: str
            ) -> None:
        ezwgd.console.log("Packup CDS, PEP and GFF3 file.")
        init_start = time.time()
        with ezwgd.console.status("Initialing..."):
            try:
                self.cds_seq = SeqIO.index(cds_file_path, "fasta")
            except FileNotFoundError:
                raise(SeqFileNotFoundError(cds_file_path))
            try:
                self.pep_seq = SeqIO.index(pep_file_path, "fasta")
            except FileNotFoundError:
                raise(SeqFileNotFoundError(pep_file_path))

            # Read and format original GFF3.
            try:
                gff = pd.read_csv(
                    gff3_file_path, sep="\t", header=None, comment="#",
                    names=["chr", "source", "type", "start", "end", "score", "strand", "phase", "attribute"],
                    dtype={"chr":"str", "type":"str", "start":"int", "end":"int", "strand":"str", "phase":"str", "attribute":"str"}
                )
            except FileNotFoundError:
                raise(GFF3FileNotFoundError(gff3_file_path))

            gff = gff.sort_values(['chr', 'start'])

            # Filtering useful features, and extract data from original "attribute" columns.
            gff = gff[gff["type"].isin(["gene", "mRNA", "transcript", "CDS"])].copy()
            gff["id"] = gff["attribute"].str.extract(r"ID=([^;]+)")
            gff["parent"] = gff["attribute"].str.extract(r"Parent=([^;]+)")

            # Save for next step.
            self.gff = gff.reset_index(drop=True)

        ezwgd.console.log(f"Initialization finished in {time.time() - init_start:3f} seconds.")

        self._build_indices()
        
    def _build_indices(self):
        ezwgd.console.log("Rebuild GFF3 tree structure.")
        # Rebuild tree structure (gene -> transcript).
        gff = self.gff
        gff_parse_start = time.time()
        with ezwgd.console.status("Rebuilding..."):
            # There is no need to rebuild from CDS.
            # -----------------------------------------------
            # transcript(tx) table (gene -> tx)
            tx = gff[gff["type"].isin(["mRNA","transcript"])].copy()
            tx = tx.rename(columns={"id":"tx_id", "parent":"gene_id"})
            tx = tx[tx["tx_id"].notna() & tx["gene_id"].notna()].loc[:,["gene_id", "tx_id"]]

            # gene -> tx dict
            self._gene_txs = tx.groupby("gene_id")["tx_id"].apply(list).to_dict()

            # -----------------------------------------------
            # gene table (for rename gene and counting genes)
            gene = gff[gff["type"] == "gene"].copy()
            gene = gene.loc[:,["id", "chr", "start", "end", "strand"]].rename(columns={"id": "gene_id"}).sort_values(["chr", "start"]).reset_index(drop=True)

            # -----------------------------------------------
            # coding gene list (as default list to extract)
            genelist = tx["gene_id"].drop_duplicates()
            gene['order'] = gene.groupby('chr').cumcount() + 1
            self.genelist = genelist.to_list()

            self.simp_gff = pd.merge(genelist, gene, on="gene_id").loc[:,['gene_id', 'chr', 'order', 'start', 'end', 'strand']]

        ezwgd.console.log(f"Rebuild finished in {time.time() - gff_parse_start:3f} seconds. "
                          f"coding gene entries: {len(self.genelist)}, all gene entries: {len(set(gene['gene_id']))}, "
                          f"transcript entries: {tx.shape[0]}.")


def fasta_extractor(
        id_list: list,
        fasta_file_path: Optional[str] = None,
        fasta_data: Optional[list[SeqRecord]] = None,
        output_path: Optional[str] = None):
    ezwgd.console.log('ezwgd.utils.parse.fasta_extractor launched.')
    
    if fasta_file_path is None and fasta_data is None:
        raise ValueError("You need input at least one of fasta_file_path or fasta_data.")
    if fasta_file_path is not None and fasta_data is not None:
        raise ValueError("You cannot input both fasta_file_path and fasta_data.")

    if fasta_file_path is not None:
        fasta_iter = SeqIO.parse(fasta_file_path, 'fasta')
        source = fasta_file_path
    else:
        fasta_iter = fasta_data
        source = 'fasta_data'

    records_dict = {record.id: record for record in fasta_iter} # type: ignore
    
    res_buffer = []
    for id_name in track(id_list, f'Extracting from {source}...'):
        try:
            res_buffer.append(records_dict[id_name])
        except KeyError:
            ezwgd.console.log(f'ID {id_name} not found in {source}, skipping.')
    if output_path is not None:
        ezwgd.console.log(f'Sequence data saved as {output_path}.')
        SeqIO.write(res_buffer, output_path, 'fasta')
    return res_buffer


def blast6reader(blast_6_result_path: str) -> pd.DataFrame:
    res = pd.read_csv(
        blast_6_result_path, sep="\t", header=None, comment="#",
        names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore",],
		dtype={"qseqid": str,
         	   "sseqid": str,
               "pident": float,
               "length": int,
               "mismatch": int,
               "gapopen": int,
               "qstart": int,
               "qend": int,
               "sstart": int,
               "send": int,
               "evalue": float,
               "bitscore": float,})
    res = res.sort_values(["qseqid", "sseqid", "bitscore"], ascending=[True, True, False])
    return res

def mcscan_res_reader(mcscan_res_path: str) -> pd.DataFrame:
    with open(mcscan_res_path, 'r') as mcscan_res, tempfile.NamedTemporaryFile('w+') as tmpf:
        bkinfo = []
        for line in mcscan_res:
            if mcscan_info_line.match(line):
                bkinfo = list(mcscan_info_line.match(line).groups()) # type: ignore
                continue
            else:
                n = line.strip().split()
                n.extend(bkinfo)
                # n has 12 columns
                #     0      1     2      3         4     5     6      7 8    9   10              11
                # gene1 order1 gene2 order2 direction bkidx score pvalue N chr1 chr2 total_direction
                # return dataframe:
                # bkidx total_direction score pvalue gene1 chr1 order1 gene2 chr2 order2 direction
                # which is
                # [5, 11, 6, 7, 0, 9, 1, 2, 10, 3, 4]
                tmpf.write('\t'.join([n[5], n[11], n[6], n[7], n[0], n[9], n[1], n[2], n[10], n[3], n[4]]) + '\n')
        tmpf.seek(0)
        return pd.read_csv(
            tmpf.name,
            sep='\t',
            header=None,
            names=['bkidx', 'total_direction', 'score', 'pvalue', 'gene1', 'chr1', 'order1', 'gene2', 'chr2', 'order2', 'direction',],
            dtype={'bkidx': int, 'total_direction': str, 'score': int, 'pvalue': float, 'gene1': str, 'chr1': str, 'order1': int,
                   'gene2': str, 'chr2': str, 'order2': int, 'direction': int}
                   )

def codeml_pairwise(rst_path: str) -> dict:
    res_name = ['N', 'S', 'dN', 'dS', 'omega', 't']

    with open(rst_path, encoding='utf-8', errors='replace') as rst:
        for line in rst.readlines():
            if pairwise_re.match(line):
                return dict(zip(res_name, pairwise_re.match(line).groups()))  # type: ignore
    raise ValueError('internal error: falled to parse result of pairwise')


def yn00_result(res_path: str):
    res_name = ['S', 'N', 't', 'kappa', 'omega', 'dN', 'dNSE', 'dS', 'dSSE']

    with open(res_path, encoding='utf-8', errors='replace') as res:
        for line in res.readlines():
            line = line.strip()
            if yn00_res_re.match(line):
                return(dict(zip(res_name, yn00_res_re.match(line).groups())))  # type: ignore
    print(line)
    raise ValueError('internal error: falled to parse result of yn00')
