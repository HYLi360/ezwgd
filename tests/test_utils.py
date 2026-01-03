"""
Testing pipeline for ezwgd module `utils`.
"""
# -------------------------------------------------------------------------------------------------------------------------
from ezwgd._custom import (
    SeqFileNotFoundError,
    GFF3FileNotFoundError,
)

from pathlib import Path
from ezwgd.utils import parse
import pytest

testing1_fna_path = "test_dataset/utils/testing1.fna"
testing1_gff3_path = "test_dataset/utils/testing1.gff3"

cds_seq0 =("ATGGTTCTCTCAAAAACTCCTTCTGATGATTCTGTACACTCCACATTTGCTTCTCGCTATGTTCGAACTT"
"CACTACCAAGGTTTGAGATGCTAGAGAAGTCTATACCAAAAGAGGCAGCATACCAAATGATTAATGATGA"
"GTTAATGCTTGATGGGAATCCAAGGTTAAATTTGGCATCATTTGTAACCACATGGATGGAACCAGAATGT"
"GATAAGCTTATGATGGCTTCAATTAACAAGAATTATGTTGACATGGATGAATACCCTGTCACCACTGAGC"
"TTCAGAATCGATGTGTAAACATGATAGCGCGTTTATTCAATGCGCCTTTGAAAGAGGAAGAAATAGGAAT"
"TGGTGTGGGGACAGTGGGGTCATCAGAGGCCATAATGTTAGCAGGGCTGGCCTTCAAGAGGAACTGGCAA"
"AACAAACGCAAAGCTGAGGGAAAGCCTTATGATAAGCCCAACATTGTCACTGGTGCTAATGTTCAGGTGT"
"GTTGGGAGAAATTTGCAAACTACTTTGAAGTGGAATTGAAACAAGTCAAGTTAAGTGAAGGGTACTATGT"
"GATGGACCCAATCAAAGCTGTGGAAATGGTAGATGACAACACTATTTGTGTTGCTGCTATTTTGGGTTCA"
"ACACTTAATGGAGAATTTGAAGATGTCAAACTCTTGAATGATCTTTTGATTGAAAAGAATAAACAAACTG"
"gATGGGACACACCTATTCATGTGGATGCAGCAAGTGGTGGATTCATTGCACCATTTATCTATCCAGAGTT"
"GGAATGGGATTTTAGGCTTCCTTTAGTGAAAAGTATTAATGTGAGTGGACACAAATATGGGCTTGTTTAT"
"GCTGGTATTGGTTGGGTTATTTGGAGAACTAAACAAGACTTGCCTCAacaactcatttttcatatcaatt"
"atCTTGGTGCTGATCAGCCTACTTTTACTCTCAATTTCTCTAAAGGTTCAAGTCAAGTCATTGCTCAATA"
"TTATCAGCTTATCCGCTTGGGCTATGAGggATATCGAAATGTAATGGAAAATTGTCGTGAAAATGCAATT"
"GTGCTAAGAAAAGGACTTGAAAAAACAGGACGTTTCAATATAATCTCCAAAGATGAAGGTATACCCTTGG"
"TGGCATTTTCCCTCAAAGACAATAGCCTCCACAACGAATTCGAGGTCTCTGAGACCCTCCGTAGGTTTGG"
"GTGGATTGTCCCAGCCTACACTATGCCAGCTGACCTGCAACATGTTACAGTGTTGCGCGTTGTGATTAGA"
"GAGGACTTCTCCCGAACCCTAGCAGATCGTCTTGTCTCTGACATCGTCAAGGTCCTCCACGAGCTCCCGA"
"ATGCCAAAAAAGTGGAGGATAATTTGATGATCAATAATGAGAAGAAAACAGAAATTGAAGTTCAAAGGGC"
"AATTGCTGAGTTTTGGAAGAAATATGTTTTAGCTAGGAAAGCATCTATTTGTTAG")

cds_seq1 = ("ATGACACATTTggaacatacatacatactcTTCTCTGCATCACACAAAAACTtatcaatcaaaaaaaaaa"
"aaatgggtgACCTTCACTCAATTGCTTACAATAGCTTCTCTTCCCTAGGTTGTGACTTCCATGGCCTAGA"
"AGTTCTCAACAATGATATGTCTCCATttatggAGAGCACTCCAccatttttgtcaaattttgatTGTGAG"
"TTGTCAAGTGGATATCTACAAGATGCTTTATTCAAATTCAACTCAAAAAGAAGgcgtttatttttgttta"
"atgatgatgatgacaaggagaattatcaaaataaagattcaattaagAATTTGTGGAGCTCCACAATTGA"
"CCAACAATTTTCTGAGGATTATGATTCTTTTAGCCAGATAACCAAGTGTGATAGCTTTTCTGgGGATCCA"
"ATGAGCAAAATGAGTGAAGAATATTCCAAAAGAACAGAGGAGGAAGCAATCtcagaaaattattattatt"
"cttctaATTCTTCACCAACATCATCATCTCATCATAATAAACAACCTCTACAATATGGAGGAGGTGAtaa"
"aaagagaagaatgtgTGGAAAAATTGTGTATCCATTTGGGCTAGTGAAGCCCGGAGGGCAAGAAGGAGAC"
"GTGACATTAAATGACATTAACGAGAGGATTCTAATGAGGCCGACCCGGCCCGTAAGGCATCCGGTTGGGG"
"ACTTTgcatgtcgttcgaccgtatgccCTGCTGGACCAGGATTATCCGGCAAGACGGTGGTTGCACTCAC"
"CAAAATTCATACTCAAGGGAGGGGGACAATCACAATTATAAGGACTCGAGGGTGA")

pep_seq0 = ("MVLSKTPSDDSVHSTFASRYVRTSLPRFEMLEKSIPKEAAYQMINDELMLDGNPRLNLASFVTTWMEPEC"
"DKLMMASINKNYVDMDEYPVTTELQNRCVNMIARLFNAPLKEEEIGIGVGTVGSSEAIMLAGLAFKRNWQ"
"NKRKAEGKPYDKPNIVTGANVQVCWEKFANYFEVELKQVKLSEGYYVMDPIKAVEMVDDNTICVAAILGS"
"TLNGEFEDVKLLNDLLIEKNKQTGWDTPIHVDAASGGFIAPFIYPELEWDFRLPLVKSINVSGHKYGLVY"
"AGIGWVIWRTKQDLPQQLIFHINYLGADQPTFTLNFSKGSSQVIAQYYQLIRLGYEGYRNVMENCRENAI"
"VLRKGLEKTGRFNIISKDEGIPLVAFSLKDNSLHNEFEVSETLRRFGWIVPAYTMPADLQHVTVLRVVIR"
"EDFSRTLADRLVSDIVKVLHELPNAKKVEDNLMINNEKKTEIEVQRAIAEFWKKYVLARKASIC")

pep_seq1 = ("MTHLEHTYILFSASHKNLSIKKKKMGDLHSIAYNSFSSLGCDFHGLEVLNNDMSPFMESTPPFLSNFDCE"
"LSSGYLQDALFKFNSKRRRLFLFNDDDDKENYQNKDSIKNLWSSTIDQQFSEDYDSFSQITKCDSFSGDP"
"MSKMSEEYSKRTEEEAISENYYYSSNSSPTSSSHHNKQPLQYGGGDKKRRMCGKIVYPFGLVKPGGQEGD"
"VTLNDINERILMRPTRPVRHPVGDFACRSTVCPAGPGLSGKTVVALTKIHTQGRGTITIIRTRG")

id0, name0, description0 = 'gene-GAD3', 'gene-GAD3', 'tx_id=rna-NM_001246898.2'
id1, name1, description1 = 'gene-LOC101263636', 'gene-LOC101263636', 'tx_id=rna-XM_004228713.5'

cds_content = f">{id0} {description0}\n{cds_seq0}\n>{id1} {description1}\n{cds_seq1}\n"
pep_content = f">{id0} {description0}\n{pep_seq0}\n>{id1} {description1}\n{pep_seq1}\n"


# -------------------------------------------------------------------------------------------------------------------------
def test_genome_general(tmp_path):
    genome = parse.Genome(
        gff3_file_path=testing1_gff3_path,
        fna_file_path=testing1_fna_path)
    
    # 1. is genome sequence exists?
    assert genome.genome_seq['seq_frag'] is not None

    # 2. is gff dataframe exists?
    assert genome.gff.shape == (17, 11)

    # 3. is get_cds_pep work correctly?
    cds, pep = genome.get_cds_pep()
    assert len(cds) == 2
    assert len(pep) == 2

    assert cds[0].seq._data == cds_seq0.encode('utf-8')
    assert cds[0].id == id0
    assert cds[0].name == name0
    assert cds[0].description == description0

    assert cds[1].seq._data == cds_seq1.encode('utf-8')
    assert cds[1].id == id1
    assert cds[1].name == name1
    assert cds[1].description == description1

    assert pep[0].seq._data == pep_seq0.encode('utf-8')
    assert pep[0].id == id0
    assert pep[0].name == name0
    assert pep[0].description == description0

    assert pep[1].seq._data == pep_seq1.encode('utf-8')
    assert pep[1].id == id1
    assert pep[1].name == name1
    assert pep[1].description == description1

    # 4. try to write file
    genome.get_cds_pep(f"{tmp_path}/cds", f"{tmp_path}/pep", warp=99999)
    assert Path(f"{tmp_path}/cds").exists()
    assert Path(f"{tmp_path}/pep").exists()
    assert Path(f"{tmp_path}/cds").read_text(encoding="utf-8") == cds_content
    assert Path(f"{tmp_path}/pep").read_text(encoding="utf-8") == pep_content


def test_genome_expertion():
    # 1. if file not exists
    with pytest.raises(SeqFileNotFoundError):
        parse.Genome("", testing1_gff3_path)

    with pytest.raises(GFF3FileNotFoundError):
        parse.Genome(testing1_fna_path, "")
    
    # 2. if index not hit
    genome = parse.Genome(
        gff3_file_path=testing1_gff3_path,
        fna_file_path=testing1_fna_path)
    
    # 2.1  Genome.get_cds(self, tx_id)
    with pytest.raises(KeyError, match=r"Transcript.+has no CDS features"):
        genome.get_cds("something")
    # compare
    cds = genome.get_cds("rna-NM_001246898.2")

    # 2.2  Genome.get_longest(self, gene_id, to_stop, strict)
    with pytest.raises(KeyError, match=r"Gene.+has no transcripts"):
        genome.get_longest("something")
    # compare
    tx_id, cds, pep = genome.get_longest("gene-GAD3")

