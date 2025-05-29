"""
Demo  radial gradient recalled echo pulse sequence

Editting parameters to obtain TR = 15ms, single slice acquisition

Test space to make changes and trial ideas
"""

import numpy as np

import pypulseq as pp


def main(plot: bool = False, write_seq: bool = False, seq_filename: str = 'gre_radial_pypulseq.seq'):
    # ======
    # SETUP
    # ======
     # FOV for SIEMENS prisma = 125 (3D), 250 (2D)
    fov = 250 # initial val =260e-3
    Nx = 64  # Define FOV and resolution
    alpha = 90 #10 = initial value # Flip angle
    slice_thickness = 1e-3 #initial val = 3e-3  # Slice thickness
    TE = 5e-3 #8e-3 = initial val  # Echo time
    TR = 7.8e-3 #20e-3 = initial val  # Repetition time
    Nr = 4 #initial val = 60  # Number of radial spokes
    N_dummy = 0  #20 = initial val # Number of dummy scans
    delta = np.pi / Nr  # Angular increment

    rf_spoiling_inc = 117  # RF spoiling increment

    # Set 3T Siemens PRISMA system limits
    system = pp.Opts(
        max_grad=139, #initial val = 28
        grad_unit='mT/m',
        max_slew=346,
        slew_unit='T/m/s',
        # Currently do not have access to these values (27/05/25)
        rf_ringdown_time=20e-6,
        rf_dead_time=100e-6,
        adc_dead_time=10e-6,
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

    ok, error_report = seq.check_timing()
    if ok:
        print('Timing check passed successfully')
    else:
        print('Timing check failed! Error listing follows:')
        print(error_report)

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


if __name__ == '__main__':
    main(plot=True, write_seq=True)
