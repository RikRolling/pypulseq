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

seq.calculate_gradient_spectrum(
    window_width=0.5,
    acoustic_resonances=mech_res,
    combine_mode='max'
     )

plt.show()
