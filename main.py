from ezwgd.frontend import ConfigCODEML, CODEML
import tempfile

with open('test_dataset/codon.aln') as codon:
    data = "".join(codon.readlines())

with tempfile.TemporaryDirectory() as tmpdir:
    codeml_res = CODEML(verbose=0).preset_dnds(tmpdir, data)
print(codeml_res)
