"""
Testing pipeline for ezwgd module `frontend`.
"""
# -------------------------------------------------------------------------------------------------------------------------
from ezwgd._custom import (
    SeqFileNotFoundError,
    GFF3FileNotFoundError,
)

from ezwgd.frontend import CODEML

def test_config():
    # frontend.CODEML
    codeml_p = CODEML(model=0, NSsites=0)
    
    codeml_p.alpha = 0.5
    assert isinstance(codeml_p.alpha, float)
    assert codeml_p.alpha == 0.5

    for v in [True, False]:
        codeml_p.cleandata = v
        assert isinstance(codeml_p.cleandata, bool)
        assert codeml_p.cleandata == v
