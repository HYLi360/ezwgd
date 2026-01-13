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

testing1_fna_path = "dataset_test/utils/testing1.fna"
testing1_gff3_path = "dataset_test/utils/testing1.gff3"

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

BLAST_RES_PATH = "dataset_test/utils/blast_res.f6"


MCSCAN_RES_PATH = "dataset_test/utils/mcscan_res.txt"

mcscan_res_target = """,bkidx,total_direction,score,pvalue,gene1,chr1,order1,gene2,chr2,order2,direction
0,12,plus,153,0.1094,Pveri1g1775,1,1775,Pvulg1g2438,1,2438,-1
1,12,plus,153,0.1094,Pveri1g1778,1,1778,Pvulg1g2440,1,2440,-1
2,12,plus,153,0.1094,Pveri1g1779,1,1779,Pvulg1g2441,1,2441,-1
3,12,plus,153,0.1094,Pveri1g1783,1,1783,Pvulg1g2443,1,2443,-1
4,12,plus,153,0.1094,Pveri1g1784,1,1784,Pvulg1g2446,1,2446,-1
5,13,plus,149,0.2704,Pveri1g597,1,597,Pvulg1g3657,1,3657,1
6,13,plus,149,0.2704,Pveri1g620,1,620,Pvulg1g3659,1,3659,-1
7,13,plus,149,0.2704,Pveri1g621,1,621,Pvulg1g3660,1,3660,-1
8,13,plus,149,0.2704,Pveri1g625,1,625,Pvulg1g3684,1,3684,1
9,13,plus,149,0.2704,Pveri1g626,1,626,Pvulg1g3689,1,3689,1
10,14,plus,137,0.1674,Pveri1g1770,1,1770,Pvulg1g2441,1,2441,-1
11,14,plus,137,0.1674,Pveri1g1772,1,1772,Pvulg1g2443,1,2443,-1
12,14,plus,137,0.1674,Pveri1g1773,1,1773,Pvulg1g2444,1,2444,-1
13,14,plus,137,0.1674,Pveri1g1777,1,1777,Pvulg1g2445,1,2445,-1
14,14,plus,137,0.1674,Pveri1g1778,1,1778,Pvulg1g2446,1,2446,-1
15,15,plus,128,0.1591,Pveri1g1700,1,1700,Pvulg1g2517,1,2517,-1
16,15,plus,128,0.1591,Pveri1g1701,1,1701,Pvulg1g2518,1,2518,1
17,15,plus,128,0.1591,Pveri1g1703,1,1703,Pvulg1g2519,1,2519,-1
18,15,plus,128,0.1591,Pveri1g1706,1,1706,Pvulg1g2522,1,2522,-1
19,15,plus,128,0.1591,Pveri1g1711,1,1711,Pvulg1g2528,1,2528,-1"""


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


def test_parse_blast6():
    blast_res = parse.blast6reader(BLAST_RES_PATH)
    # 1. is result reads?
    assert blast_res.shape == (96777, 12)

    # 2. 


def test_parse_mcscan_res():
    mcscan_res = parse.mcscan_res_reader(MCSCAN_RES_PATH)
    # 
    assert mcscan_res.shape == (20, 11)
