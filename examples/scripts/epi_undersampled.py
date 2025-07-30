"""
Demo undersampled low-performance EPI sequence without ramp-sampling.

Editting parameters to obtain TR = 15ms, single slice acquisition

Test space to make changes and trial ideas
"""

import numpy as np

import pypulseq as pp


def main(plot: bool = False, write_seq: bool = False, seq_filename: str = 'epi_undersampled_lowsysvals_longtime2.seq'):
    # ======
    # SETUP
    # ======
    # Define FOV and resolution
    # FOV for SIEMENS prisma = 125 (3D), 250 (2D)
    fov = 250e-3 #125e-3 = 3D Siemens value #220e-3 = initial val
    Nx = 64
    Ny = 64
    slice_thickness = 3e-3  # Slice thickness
    n_slices = 1
    #added undersampling factor
    undersampling_factor = 0.5

    # Set 3T Siemens PRISMA system limits
    system = pp.Opts(
        max_grad=75,
        grad_unit='mT/m',
        max_slew=280,
        slew_unit='T/m/s',
        # Currently do not have access to these values (23/05/25)
        rf_ringdown_time=30e-4,
        rf_dead_time=100e-4,
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
    #added undersampling factor
    readout_time = undersampling_factor*Nx * dwell_time
    flat_time = np.ceil(readout_time * 1e5) * 1e-5  # round-up to the gradient raster
    #added
    gx = pp.make_trapezoid(
        channel='x',
        system=system,
        amplitude=k_width / readout_time,
        flat_time=flat_time,
    )
    adc = pp.make_adc(
        num_samples=Nx*undersampling_factor,
        duration=readout_time,
        delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2,
    )

    # Pre-phasing gradients
    pre_time = 8e-4
    gx_pre = pp.make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
    gz_reph = pp.make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
    gy_pre = pp.make_trapezoid(channel='y', system=system, area=-Ny / 2 * delta_k, duration=pre_time)

    # Phase blip in the shortest possible time
    dur = np.ceil(2 * np.sqrt(delta_k / system.max_slew) / 10e-6) * 10e-6
    dur_u = np.ceil(2*2 * np.sqrt(delta_k / system.max_slew) / 10e-6) * 10e-6
    gy = pp.make_trapezoid(channel='y', system=system, area=delta_k, duration=dur)
    gy_amplitude_undersample = 2*gy.amplitude
    gy_undersample = pp.make_trapezoid(channel='y', system=system, amplitude=gy_amplitude_undersample, duration=dur_u)
    N_u = int(Ny/2)
    

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
            seq.add_block(gy_undersample)  # Phase blip
            gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

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
        seq.plot()  # Plot sequence waveforms

    # =========
    # WRITE .SEQ
    # =========
    if write_seq:
        seq.write(seq_filename)

    return seq


if __name__ == '__main__':
    main(plot=True, write_seq=True)