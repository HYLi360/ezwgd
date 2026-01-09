from ezwgd.frontend import ConfigCODEML, ConfigCODEMLAdvanced



user_cfg = ConfigCODEML(runmode='usertree', seqtype='codon', model='free-ratios', NSsites='one-omega', noisy='none')
adv_cfg = ConfigCODEMLAdvanced.escape(user_cfg)
ctl_text = str(adv_cfg)
print(ctl_text)
