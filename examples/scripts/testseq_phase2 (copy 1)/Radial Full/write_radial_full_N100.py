"""
Demo golden angle radial gradient recalled echo pulse sequence

Editting parameters to obtain TR <= 15ms, single slice acquisition

Test space to make changes and trial ideas

Half spoke, delta ~ 137.51 degrees, Gx area = Nx/2 * a

Full Spoke, delta ~ 111.25 degrees, Gx area = Nx * a

where a is some other variables.
"""

import numpy as np

import pypulseq as pp

import pandas as pd

import matplotlib.pyplot as plt

from pypulseq.SAR.SAR_calc import _SAR_from_seq as SAR

from pypulseq.SAR.SAR_calc import _load_Q

from pypulseq.utils.siemens import readasc as readasc

from pypulseq.utils.siemens import asc_to_hw as asc_to_hw


def main(plot: bool = False, write_seq: bool = False, pns_check: bool = False, test_report: bool = False, sar: bool = False , acoustic_check: bool = False ,k_space: bool = False, seq_filename: str = 'gre_radial_full_N100.seq'):
    # ======
    # SETUP
    # ======
     # FOV for SIEMENS prisma = 125 (3D), 250 (2D)
    fov = 250e-3 # initial val =260e-3
    Nx = 64  # Define FOV and resolution
    alpha = 90 #10 = initial value # Flip angle
    slice_thickness = 3e-3 #initial val = 3e-3  # Slice thickness
    TE = 5e-3  #7.5 #8e-3 = initial val  # Echo time
    TR = 10e-3 #15 #20e-3 = initial val  # Repetition time
    Nr = 100 #initial val = 60  # Number of radial spokes
    N_dummy = 0  #20 = initial val # Number of dummy scans
    delta = 111.25*(np.pi/360) # Angular increment

    rf_spoiling_inc = 117  # RF spoiling increment

    # Set 3T Siemens PRISMA system limits
    system = pp.Opts(
        max_grad=50, #initial val = 28
        grad_unit='mT/m',
        max_slew=130,#120
        slew_unit='T/m/s',
        # Currently do not have access to these values (27/05/25)
        rf_ringdown_time=20e-6,
        rf_dead_time=100e-6,
        adc_dead_time=10e-6
    )

    seq = pp.Sequence(system)  # Create a new sequence object

    # ======
    # CREATE EVENTS
    # ======
    # Create alpha-degree slice selection pulse and gradient
    rf, gz, _ = pp.make_sinc_pulse(
        apodization=0.5,
        duration=4e-3,
        flip_angle=alpha * np.pi / 180,
        slice_thickness=slice_thickness,
        system=system,
        time_bw_product=4,
        return_gz=True,
        delay=system.rf_dead_time,
    )

    # Define other gradients and ADC events
    #Trial Nx/2 for half spoke
    deltak = 1 / fov
    gx = pp.make_trapezoid(channel='x', flat_area=Nx * deltak, flat_time=6.4e-3 / 5, system=system)
    adc = pp.make_adc(num_samples=Nx, duration=gx.flat_time, delay=gx.rise_time, system=system)
    gx_pre = pp.make_trapezoid(channel='x', area=-gx.area / 2 - deltak / 2, duration=2e-3, system=system)
    gz_reph = pp.make_trapezoid(channel='z', area=-gz.area / 2, duration=2e-3, system=system)
    # Gradient spoiling
    gx_spoil = pp.make_trapezoid(channel='x', area=0.5 * Nx * deltak, system=system)
    gz_spoil = pp.make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

    # Calculate timing
    delay_TE = (
        np.ceil(
            (TE - pp.calc_duration(gx_pre) - gz.fall_time - gz.flat_time / 2 - pp.calc_duration(gx) / 2)
            / seq.grad_raster_time
        )
        * seq.grad_raster_time
    )
    delay_TR = (
        np.ceil(
            (TR - pp.calc_duration(gx_pre) - pp.calc_duration(gz) - pp.calc_duration(gx) - delay_TE)
            / seq.grad_raster_time
        )
        * seq.grad_raster_time
    )
    assert np.all(delay_TR) > pp.calc_duration(gx_spoil, gz_spoil)
    rf_phase = 0
    rf_inc = 0

    # ======
    # CONSTRUCT SEQUENCE
    # ======
    for i in range(-N_dummy, Nr + 1):
        rf.phase_offset = rf_phase / 180 * np.pi
        adc.phase_offset = rf_phase / 180 * np.pi
        rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
        rf_phase = divmod(rf_inc + rf_phase, 360.0)[1]

        seq.add_block(rf, gz)
        #read through k-space defined here
        phi = delta * (i - 1)
        seq.add_block(*pp.rotate(gx_pre, gz_reph, angle=phi, axis='z'))
        seq.add_block(pp.make_delay(delay_TE))
        if i > 0:
            seq.add_block(*pp.rotate(gx, adc, angle=phi, axis='z'))
        else:
            seq.add_block(*pp.rotate(gx, angle=phi, axis='z'))
        seq.add_block(*pp.rotate(gx_spoil, gz_spoil, pp.make_delay(delay_TR), angle=phi, axis='z'))
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
        #Using combined asc files
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
        with open('test_report_radialfull_N100.txt', 'w') as file:
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
    # ========
    # SAR CHECKER
    # ========

    # USE OF SAR IS INCORRECT!!! wE REQUIRE Q MATRIX
    if sar:

        Qtmf, Qhmf = _load_Q()
        sar_values = SAR(seq, Qtmf, Qhmf)
        sar_values_array = np.column_stack((sar_values[0], sar_values[1], sar_values[2]))

        headers = ["Body mass SAR", "Head mass SAR", "time"]
        sar_values_table = pd.DataFrame(sar_values_array, columns=headers)
        sar_values_table.to_csv('SAR_fullrad_N100.csv', index=False)

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
        seq.plot()

    # =========
    # WRITE .SEQ
    # =========
    if write_seq:
        seq.set_definition(key='FOV', value=[fov, fov, slice_thickness])
        seq.set_definition(key='Name', value='gre_rad')
        seq.write(seq_filename)

    return seq

    #plt.close()

if __name__ == '__main__':
    main(plot=True, write_seq=True, pns_check=True, test_report=True, sar=True, acoustic_check=True, k_space=True)
