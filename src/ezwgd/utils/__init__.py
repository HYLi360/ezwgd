import re

pairwise_re = re.compile(r"^\s*\d+\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")
"""N S dN dS dN/dS(omega) t"""

mcscan_info_line = re.compile(r"^#\sAlignment\s(\d+):\sscore=(\d+)\spvalue=([\d.]+)\sN=(\d+)\s(\d+)&(\d+)\s(\w+)")
"""bkidx, score, pvalue, N, chr1, chr2, direction"""

yn00_res_re = re.compile(r"^2\s+1\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\+\-\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\+\-\s+\-?([\d\.]+)")
"""S N t kappa omega dN dNSE dS dSSE"""
