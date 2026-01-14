"""
A collinearity detector based on "improved collinearity" module of WGDI.

"""


from ezwgd.utils.parse import Genome, GenomeLight

import pandas as pd


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
    def chunking(self):
        """Chunk the blast6 result to prepare parallel processes."""
        blast_res = self.blast_res.loc[:,["qseqid", "sseqid", "evalue", "bitscore",]]
