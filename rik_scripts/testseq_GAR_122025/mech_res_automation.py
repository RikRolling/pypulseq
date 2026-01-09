"""
Shifts peaks outside of mechanical resonance bandwidth 
for a given input sequence.

Returns new sequence, test report with sequence parameters,
PNS report, mech res report and new gradient/adc diagram.

Author: Rik Khot
Github: @RikRolling

Work is part of 'The Beat Goes On' EPSRC phd studentship, grant no. 524682/1. (2025-2029)
"""

import numpy as np

import pypulseq as pp

import pandas as pd

import matplotlib.pyplot as plt

from pypulseq.SAR.SAR_calc import _SAR_from_seq as SAR

from pypulseq.SAR.SAR_calc import _load_Q

from pypulseq.utils.siemens import readasc as readasc

from pypulseq.utils.siemens import asc_to_hw as asc_to_hw


# Create sequence object
seq = pp.Sequence()

# User enters file path to sequence and scanner specific asc file
seq.read('/cubric/data/c24073803/pypulseq_repo/pypulseq/rik_scripts/testseq_GAR_092025/GAR_N100_TR40e-4_cons/GAR_N100_TR40e-4_cons.seq')
asc, extra = readasc.readasc('/cubric/data/c24073803/pypulseq_repo/pypulseq/rik_scripts/combined_copy.asc')


mech_res = asc_to_hw.asc_to_acoustic_resonances(asc)
print('mechanical resonances:', mech_res)

mech_res_1 = mech_res[0]
mech_res_2 = mech_res[1]

f1_res = mech_res_1['frequency']
f1_bwth = mech_res_1['bandwidth']

f2_res = mech_res_2['frequency']
f2_bwth = mech_res_2['bandwidth']

print('mech res 1',f1_res, '+/-', f1_bwth/2)
print('mech res 2',f2_res, '+/-', f2_bwth/2)

spectogram,spectogram_rss, f, t = seq.calculate_gradient_spectrum(
    window_width=0.5,
    acoustic_resonances=mech_res,
    combine_mode='max',
    plot=True
     )

plt.show()

# Define Mask
lb1 = f1_res - f1_bwth/2
ub1 = f1_res + f1_bwth/2
lb2 = f2_res - f2_bwth/2
ub2 = f2_res + f2_bwth/2


#if lb1 <= f <= ub1 or lb2 <= f <= ub2:

