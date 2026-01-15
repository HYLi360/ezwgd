"""
A collinearity detector based on "improved collinearity" module of WGDI.

"""


from ezwgd.utils.parse import Genome, GenomeLight

import pandas as pd
import numpy as np

class PowerColl:
    """
    Run collineary block calculate per chromosome pair.\n
    Input: \n
    Total gene sequence and their chromosome name.\n
    `gene_sequence = ['gene1', 'gene2', ...]`\n
    `chr_name = ['chr1', 'chr1', ...]`\n
    BLASTp results of gene pair (from BLAST++, diamond and so on, by blast6reader.)\n
    """
    def __init__(
            self,
            blast_res: pd.DataFrame,
            genome1 : Genome | GenomeLight,
            genome2 : Genome | GenomeLight,
            ) -> None:
        self.blast_res = blast_res
        self.genome1 = genome1
        self.genome2 = genome2
        self.simp_gff1 = genome1.simp_gff
        self.simp_gff2 = genome2.simp_gff
    def chunking(self):
        """Chunk the blast6 result to prepare parallel processes."""
        blast_res = self.blast_res.loc[:,["qseqid", "sseqid", "evalue", "bitscore",]]
        blast_res.columns = ['gene_id1', 'gene_id2', 'evalue', 'bitscore']
        gene2chr1, gene2chr2 = dict(zip(self.simp_gff1['gene_id'], self.simp_gff1['chr'])), dict(zip(self.simp_gff2['gene_id'], self.simp_gff2['chr']))
        gene2order1, gene2order2 = dict(zip(self.simp_gff1['gene_id'], self.simp_gff1['order'])), dict(zip(self.simp_gff2['gene_id'], self.simp_gff2['order']))
        blast_res['chr_id1'] = blast_res['gene_id1'].map(gene2chr1)
        blast_res['chr_id2'] = blast_res['gene_id2'].map(gene2chr2)
        blast_res['order1'] = blast_res['gene_id1'].map(gene2order1)
        blast_res['order2'] = blast_res['gene_id2'].map(gene2order2)

        blast_res = blast_res.loc[:,['order1', 'order2', 'evalue', 'bitscore']]

        blast_res = blast_res.sort_values(['order1', 'order2']).reset_index(drop=True).loc[:,['order1', 'order2']]
        blast_res.columns = ['loc1', 'loc2']
        blast_res['grading'] = 50

        blast_res.to_csv("res.tsv", sep='\t')