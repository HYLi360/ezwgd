"""
Testing pipeline for ezwgd module `frontend.PAML`.
"""
# -------------------------------------------------------------------------------------------------------------------------
from dataclasses import fields

from ezwgd.frontend.PAML import ConfigCODEML, CODEML

def test_config_default_setting():
    config = ConfigCODEML()
    config_advanced = CODEML()
    param_list = ['seqfile', 'outfile', 'noisy', 'verbose', 'runmode', 'seqtype', 'CodonFreq', 'ndata', 'clock', 'model',
                  'NSsites', 'icode', 'Mgene', 'fix_kappa', 'kappa', 'fix_omega', 'omega', 'fix_alpha', 'alpha',
                  'Malpha', 'fix_rho', 'rho', 'ncatG', 'getSE', 'RateAncestor', 'Small_Diff', 'cleandata', 'fix_blength', 'method']
    param_default_semantic = {'seqfile': '',
                              'outfile': '',
                              'noisy': 'max',
                              'verbose': False,
                              'runmode': 'usertree',
                              'seqtype': 'aa',
                              'CodonFreq': 'F3X4',
                              'ndata': 10,
                              'clock': 'no',
                              'model': 'more-per-branch',
                              'NSsites': 'one-omega',
                              'icode': 'Standard',
                              'Mgene': 'rates',
                              'fix_kappa': False,
                              'kappa': 2,
                              'fix_omega': False,
                              'omega': 0.4,
                              'fix_alpha': True,
                              'alpha': 0,
                              'Malpha': False,
                              'ncatG': 3,
                              'fix_rho': True,
                              'rho': 0,
                              'getSE': False,
                              'RateAncestor': 0,
                              'Small_Diff': 5e-7,
                              'cleandata': False,
                              'fix_blength': 'ignore',
                              'method': 'simultaneous'
                              }
    param_default_numeric = {'seqfile': '',
                              'outfile': '',
                              'noisy': 9,
                              'verbose': 0,
                              'runmode': 0,
                              'seqtype': 2,
                              'CodonFreq': 2,
                              'ndata': 10,
                              'clock': 0,
                              'model': 2,
                              'NSsites': 0,
                              'icode': 0,
                              'Mgene': 0,
                              'fix_kappa': 0,
                              'kappa': 2,
                              'fix_omega': 0,
                              'omega': 0.4,
                              'fix_alpha': 1,
                              'alpha': 0,
                              'Malpha': 0,
                              'ncatG': 3,
                              'fix_rho': 1,
                              'rho': 0,
                              'getSE': 0,
                              'RateAncestor': 0,
                              'Small_Diff': 5e-7,
                              'cleandata': 0,
                              'fix_blength': 0,
                              'method': 0
                              }
    # Verify that the default settings for `ConfigCODEML` and `CODEML` are correct and consistent.
    # 1. parameters name check
    # `ConfigCODEML` and `CODEML` should have 29 parameters (exclude `treefile`, which is optionable conditionally)
    assert set([param.name for param in fields(config)]) == set(param_list)
    assert set([param.name for param in fields(config_advanced)]) == set(param_list)
    # Ensure that the keys in the dictionary used for verification are also consistent
    assert set(param_default_semantic.keys()) == set(param_list)
    assert set(param_default_numeric.keys()) == set(param_list)

    # 2. default parameters check
    for param in param_list:
        assert param_default_semantic[param] == getattr(config, param)
        assert param_default_numeric[param] == getattr(config_advanced, param)
    
    # 3. escape functions
    config_advanced = CODEML.escape(config)
    assert set([param.name for param in fields(config_advanced)]) == set(param_list)
    for param in param_list:
        assert param_default_numeric[param] == getattr(config_advanced, param)

    # 4. other parameters

    # 5. actual functions
