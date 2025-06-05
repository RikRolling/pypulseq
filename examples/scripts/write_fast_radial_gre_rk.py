"""
Translation of MATLAB pulseq tutorial.
Original MATLAB author: Maxim Zaitsev
MATLAB to Python translation author: Rik Khot

Fast Radial Gradient Echo Pulse Sequence

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
    slice_thickness = 1e-3 #initial val = 3e-3  

    #Use above parameters to find minimum TE and TE
   
    Nr = 60 #initial val = 60  # Number of radial spokes
    N_dummy = 20  #20 = initial val # Number of dummy scans
    delta = np.pi / Nr  # Angular increment

    # readout parameters

    ro_dur = 1200e-6 #readout duration
    ro_os = 2 # oversampling factor
    ro_spoil = 0.5
    rf_dur = 600e-6
    sl_spoil = 2


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
    rf, gz, gz_reph = pp.make_sinc_pulse(
        apodization=0.5,
        duration=rf_dur,
        flip_angle=alpha * np.pi / 180,
        slice_thickness=slice_thickness,
        system=system,
        time_bw_product=2,
        return_gz=True,
    )
    gz_reph.delay = pp.calc_duration(gz)
    gz_comb = pp.add_gradients(grads=[gz, gz_reph],system=system)

    # Define other gradients and ADC events
    deltak = 1 / fov
    gx = pp.make_trapezoid(channel='x', amplitude=Nx * deltak/ ro_dur, flat_time=np.ceil(ro_dur/system.grad_raster_time)*system.grad_raster_time, system=system)
    adc = pp.make_adc(num_samples=Nx, duration=ro_dur, delay=gx.rise_time, system=system)
    gx_pre = pp.make_trapezoid(channel='x', area=-gx.amplitude*(ro_dur/Nx/ro_os*(Nx*ro_os/2-0.5)+0.5*gx.rise_time),system=system) #0.5 used to account for Siemens sampling at centre of dwell periods
    gx_pre,_,_ = pp.align(right=[gx_pre,gz_comb], left=pp.make_delay(pp.calc_duration(rf)+pp.calc_duration(gx_pre)))
   
    # Gradient spoiling

    if sl_spoil > 0:
        sp_area_needed=sl_spoil/slice_thickness-gz.area/2
        gz_spoil=pp.make_trapezoid(channel='z', area=sp_area_needed, system=system, delay=gx.rise_time + gx.flat_time)
    else:
        gz_spoil=[]

    if ro_spoil > 0:
        ro_add_time = np.ceil(((gx.area/Nx*(Nx/2+1)*ro_spoil)/gx.amplitude)/system.grad_raster_time)*system.grad_raster_time
        gx.flat_time = gx.flat_time + ro_add_time

    #We accept what time is achievable, time is not calculated

    TR = 0
    TE = 0

    rf_phase = 0
    rf_inc = 0

    duration, num_blocks, event_count = iter(seq.duration())
    #seq_duration_tuple = tuple(seq.duration)

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
        seq.add_block(*pp.rotate(rf, gz_comb, gx_pre, angle=phi, axis='z'))

        if TE<=0:
            

            TE = duration + adc.delay + adc.dwell*(adc.num_samples/2+0.5)
       
        if i > 0:
            seq.add_block(*pp.rotate(gx, adc, gz_spoil, angle=phi, axis='z'))

        else:
            seq.add_block(*pp.rotate(gx, gz_spoil, angle=phi, axis='z'))
        
        if TR<=0:
            TR=duration
      

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
