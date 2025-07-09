"""
Demo low-performance EPI sequence without ramp-sampling.

RIK KHOT

Fully Sampled EPI, TR ~ 40ms
"""

import numpy as np

import pandas as pd

import pypulseq as pp

import matplotlib.pyplot as plt

import scipy as sp

from scipy.optimize import minimize, LinearConstraint

from pypulseq.SAR.SAR_calc import _SAR_from_seq as SAR

from pypulseq.SAR.SAR_calc import _load_Q

from pypulseq.utils.siemens import readasc as readasc

from pypulseq.utils.siemens import asc_to_hw as asc_to_hw

from pypulseq import calc_rf_bandwidth

from pypulseq import calc_rf_center


def main(plot: bool = False, write_seq: bool = False, pns_check: bool = False, test_report: bool = False, sar: bool = False , acoustic_check: bool = False ,k_space: bool = False, seq_filename: str = 'epi_fullsample.seq'):
    # ======
    # SETUP
    # ======
    # Define FOV and resolution
    fov = 250e-3
    Nx = 64
    Ny = 64
    slice_thickness = 3e-3  # Slice thickness
    n_slices = 1

    # Set system limits
    system = pp.Opts(
        max_grad=80,
        grad_unit='mT/m',
        max_slew=150,
        slew_unit='T/m/s',
        rf_ringdown_time=30e-6,
        rf_dead_time=100e-6,

    )

    seq = pp.Sequence(system)  # Create a new sequence object

    # ======
    # CREATE EVENTS
    # ======
    # Create 90 degree slice selection pulse and gradient
    rf, gz, _ = pp.make_sinc_pulse(
        flip_angle=np.pi / 2,
        system=system,
        duration=3e-3,
        slice_thickness=slice_thickness,
        apodization=0.5,
        time_bw_product=4,
        return_gz=True,
        delay=system.rf_dead_time,
    )

    # Define other gradients and ADC events
    delta_k = 1 / fov
    k_width = Nx * delta_k
    dwell_time = 4e-6
    readout_time = Nx * dwell_time
    flat_time = np.ceil(readout_time * 1e5) * 1e-5  # round-up to the gradient raster
    #print(flat_time)
    gx = pp.make_trapezoid(
        channel='x',
        system=system,
        amplitude=k_width / readout_time,
        flat_time=flat_time,
    )
    adc = pp.make_adc(
        num_samples=Nx,
        duration=readout_time,
        delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2,
    )


    # Pre-phasing gradients
    pre_time = 8e-4
    gx_pre = pp.make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
    gz_reph = pp.make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
    gy_pre = pp.make_trapezoid(channel='y', system=system, area=-Ny / 2 * delta_k, duration=pre_time)

    # Phase blip in the shortest possible time
    dur =np.ceil(2 * np.sqrt(delta_k /system.max_slew) / 10e-6) * 10e-6
    gy = pp.make_trapezoid(channel='y', system=system, area=delta_k, duration=dur)
    # ======
    # CONSTRUCT SEQUENCE
    # ======
    # Define sequence blocks
    for s in range(n_slices):
        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf, gz)
        seq.add_block(gx_pre, gy_pre, gz_reph)
        for _ in range(Ny):
            seq.add_block(gx, adc)  # Read one line of k-space
            seq.add_block(gy)  # Phase blip
            gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient


    # ======
    # Timing Check
    # ======
    ok, error_report = seq.check_timing()
    if ok:
        print('Timing check passed successfully')


    else:
        print('Timing check failed! Error listing follows:')
        print(error_report)

    # ======
    # PNS Checker
    # ======
    if pns_check:
        #Combine asc files
        a, b, c, d = seq.calculate_pns('combined_copy.asc')
        if a == True:
            print('PNS check passed')
        if a == False:
            print('PNS check failed')

    # ======
    # Test Report
    # ======
    if test_report:
        #user to change text name based on read-out trajectory
        with open('test_report_epi.txt', 'w') as file:
            file.write(seq.test_report())

    # ======
    # Accoustic Frequency Checker
    # ======
    if acoustic_check:
        asc, extra = readasc.readasc('combined_copy.asc')
        list = asc_to_hw.asc_to_acoustic_resonances(asc)
        seq.calculate_gradient_spectrum(
            acoustic_resonances=list,
            plot=True
        )


    # ======
    # SAR
    # ======
    if sar:

        Qtmf, Qhmf = _load_Q()
        #print(Qtmf)
        sar_values = SAR(seq, Qtmf, Qhmf)
        sar_values_array = np.column_stack((sar_values[0], sar_values[1], sar_values[2]))

        headers = ["Body mass SAR", "Head mass SAR", "time"]
        sar_values_table = pd.DataFrame(sar_values_array, columns=headers)
        sar_values_table.to_csv('SAR.csv', index=False)

        #SAR checker - print statement will only been shown if SAR is violated for either head or body
        violation_1 = False
        violation_2 = False
        for i in sar_values_table.iloc[:, 1]:
            if (i >= 3.2):
                print("SAR head value NOT acceptable")
                violation_1 = True
                break
        #Validation for full body mass SAR
        for j in sar_values_table.iloc[:, 0]:
            if (j > 2):
                print("SAR Body value NOT acceptable")
                violation_2 = False
                break
    # ======
    # K-SPACE VISUALIZATION
    # ======
    if k_space:

        k_adc, k_all, t_excite, t_refocus, t_adc = seq.calculate_kspace()
        plt.figure()
        plt.title('Estimated K-Space Trajectory')
        plt.plot(k_adc[0], k_adc[1])
        plt.show()






    # ======
    # VISUALIZATION
    # ======
    if plot:

        seq.plot()  # Plot sequence waveforms

    # =========
    # WRITE .SEQ
    # =========
    if write_seq:
        seq.write(seq_filename)

    return seq


if __name__ == '__main__':
    main(plot=True, write_seq=True, pns_check=True, test_report=True, sar=True, acoustic_check=True, k_space=True)