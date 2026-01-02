"""
Testing pipeline for ezwgd module `utils`.
"""

testing_gff3_path = "../../test_dataset/testing.gff3"
testing_fna_path = "../../test_dataset/testing.fna"

from ezwgd.utils import parse

def test1_genome_parse():
    genome = parse.Genome(
        gff3_file_path=testing_gff3_path,
        fna_file_path=testing_fna_path)
    
    # 1. is genome sequence exists?
    assert genome.genome_seq['seq_frag'] is not None

    # 2. 