
"""
Testing pipeline for ezwgd module `frontend`.
"""
# -------------------------------------------------------------------------------------------------------------------------

from ezwgd.frontend import ConfigCODEML


def test_config():
    # frontend.ConfigCODEML
    codeml_p = ConfigCODEML()

    codeml_p.alpha = 0.5
    assert isinstance(codeml_p.alpha, float)
    assert codeml_p.alpha == 0.5

    for v in [True, False]:
        codeml_p.cleandata = v
        assert isinstance(codeml_p.cleandata, bool)
        assert codeml_p.cleandata == v


def test_gffread_cds_pep_extract():
    assert False


def test_muscle5():
    assert False


def test_mafft():
    assert False


def test_diamond():
    assert False


def test_config_codeml():
    assert False


def test_config_yn00():
    assert False


def test_chi2():
    assert False
