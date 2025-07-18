#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple EPI Sequence with spoiled gradients

System limits based on SIEMENS 3T PRISMA

TR = 46ms
Nx, Ny = 32
TE = 37e-3
pre-time = 5e-4
rf180 duration = 100e-6
rf_dead_time = 100e-5

@author: Rik Khot
"""

import math

import numpy as np

import pypulseq as pp


def main(plot: bool = False, write_seq: bool = False, seq_filename: str = 'epi_se_46ms_pypulseq.seq'):
    # ======
    # SETUP
    # ======
    # FOV for SIEMENS prisma = 125 (3D), 250 (2D)
    fov = 250e-3  # Define FOV and resolution
    Nx = 32
    Ny = 32

    # Set system limits
    # Set 3T Siemens PRISMA system limits
    system = pp.Opts(
        max_grad=139,
        grad_unit='mT/m',
        max_slew=346,
        slew_unit='T/m/s',
        # Currently do not have access to these values (23/05/25)
        rf_ringdown_time=20e-6,
        rf_dead_time=100e-5, #100e-6 - RF dead time
        adc_dead_time=20e-6,
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
        slice_thickness=3e-3,
        apodization=0.5,
        time_bw_product=4,
        return_gz=True,
        delay=system.rf_dead_time,
    )

    # Define other gradients and ADC events
    delta_k = 1 / fov
    k_width = Nx * delta_k
    readout_time = 3.2e-4
    gx = pp.make_trapezoid(channel='x', system=system, flat_area=k_width, flat_time=readout_time)
    adc = pp.make_adc(num_samples=Nx, system=system, duration=gx.flat_time, delay=gx.rise_time)

    # Pre-phasing gradients
    pre_time = 5e-4 #8e-4 inital value
    gz_reph = pp.make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
    # Do not need minus for in-plane prephasers because of the spin-echo (position reflection in k-space)
    gx_pre = pp.make_trapezoid(channel='x', system=system, area=gx.area / 2 - delta_k / 2, duration=pre_time)
    gy_pre = pp.make_trapezoid(channel='y', system=system, area=Ny / 2 * delta_k, duration=pre_time)

    # Phase blip in shortest possible time
    dur = math.ceil(2 * math.sqrt(delta_k / system.max_slew) / 10e-6) * 10e-6
    gy = pp.make_trapezoid(channel='y', system=system, area=delta_k, duration=dur)

    # Refocusing pulse with spoiling gradients
    rf180 = pp.make_block_pulse(
        flip_angle=np.pi,
        delay=system.rf_dead_time,
        system=system,
        duration=100e-6, #500e-6 - original value
        use='refocusing',
    )
    gz_spoil = pp.make_trapezoid(channel='z', system=system, area=gz.area * 2, duration=3 * pre_time)

    # Calculate delay time
    TE = 37e-3 #original value 60e-3
    duration_to_center = (Nx / 2 + 0.5) * pp.calc_duration(gx) + Ny / 2 * pp.calc_duration(gy)
    rf_center_incl_delay = rf.delay + pp.calc_rf_center(rf)[0]
    rf180_center_incl_delay = rf180.delay + pp.calc_rf_center(rf180)[0]
    delay_TE1 = (
        TE / 2
        - pp.calc_duration(gz)
        + rf_center_incl_delay
        - pre_time
        - pp.calc_duration(gz_spoil)
        - rf180_center_incl_delay
    )
    delay_TE2 = (
        TE / 2 - pp.calc_duration(rf180) + rf180_center_incl_delay - pp.calc_duration(gz_spoil) - duration_to_center
    )

    # ======
    # CONSTRUCT SEQUENCE
    # ======
    # Define sequence blocks
    seq.add_block(rf, gz)
    seq.add_block(gx_pre, gy_pre, gz_reph)
    seq.add_block(pp.make_delay(delay_TE1))
    seq.add_block(gz_spoil)
    seq.add_block(rf180)
    seq.add_block(gz_spoil)
    seq.add_block(pp.make_delay(delay_TE2))
    for _ in range(Ny):
        seq.add_block(gx, adc)  # Read one line of k-space
        seq.add_block(gy)  # Phase blip
        gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient
    seq.add_block(pp.make_delay(1e-4))

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
        seq.write(seq_filename)

    return seq


if __name__ == '__main__':
    main(plot=True, write_seq=True)


