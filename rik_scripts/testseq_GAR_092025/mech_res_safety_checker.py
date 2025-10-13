#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 15:18:33 2025

@author: ritikakhot
"""

import numpy as np

import pypulseq as pp

import pandas as pd

import matplotlib.pyplot as plt

from pypulseq.SAR.SAR_calc import _SAR_from_seq as SAR

from pypulseq.SAR.SAR_calc import _load_Q

from pypulseq.utils.siemens import readasc as readasc

from pypulseq.utils.siemens import asc_to_hw as asc_to_hw


seq = pp.Sequence()
seq.read('/Users/ritkakhot/GitHub/rik_scripts/testseq_GAR_092025/GAR_N100_TR40e-4/GAR_N100_TR40e-4.seq')

asc, extra = readasc.readasc('combined_copy.asc')
list = asc_to_hw.asc_to_acoustic_resonances(asc)
spectogram, spectogram_sos, f, t = seq.calculate_gradient_spectrum(
    acoustic_resonances=list,
    combine_mode='max',
    plot=True
     )

seq.calculate_gradient_spectrum(
acoustic_resonances=list,
combine_mode='max',
plot=True
 )

print("Spectrogram shape:", spectogram_sos.shape)
print("Frequency array shape:", f.shape)


# --- Step 2: Filter by resonance bands ---
bands = [(490, 690), (940, 1360)]
mask = np.zeros_like(f, dtype=bool)
for low, high in bands:
    mask |= (f >= low) & (f <= high)

# --- Step 3: Compute total and filtered integrals ---
# Assuming spectrogram_sos is 1D (same length as f)
if spectogram_sos.ndim == 1:
    total_integral = np.trapz(spectogram_sos, f)
    filtered_integral = np.trapz(spectogram_sos[mask], f[mask])
else:
# If spectrogram_sos is 2D (e.g., time × frequency)
    total_integral = np.trapz(np.trapz(spectogram_sos, f, axis=1), dx=1)
    filtered_integral = np.trapz(np.trapz(spectogram_sos[:, mask], f[mask], axis=1), dx=1)

# --- Step 4: Compute ratio ---
percent_filtered = (filtered_integral / total_integral) * 100

print(f"Total spectrogram integral: {total_integral:.4f}")
print(f"Filtered spectrogram integral: {filtered_integral:.4f}")
print(f"Mechanical resonance ratio: {percent_filtered:.2f}%")

# Optional: store filtered subset if you want
filtered_data = np.vstack([f[mask], spectogram_sos[mask]]).T
