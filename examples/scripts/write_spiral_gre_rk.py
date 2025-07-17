"""
Converting writeSpiral.m from matlab to python

Original author: Maxim Zaitsev

Author for .py: Rik Khot
"""

import numpy as np

import pypulseq as pp

import pandas as pd

import matplotlib as plt

from pypulseq.SAR.SAR_calc import _SAR_from_seq as SAR

<<<<<<< HEAD
from pypulseq.SAR.SAR_calc import _load_Q 
 
def main(plot: bool = False, write_seq: bool = False, pns_check: bool = False, test_report: bool = False, sar: bool = False , acoustic_check: bool = False ,k_space: bool = False, seq_filename: str = 'gre_spiral.seq'):
=======
from pypulseq.SAR.SAR_calc import _load_Q

def main(plot: bool = False, write_seq: bool = False, sar: bool = False , seq_filename: str = 'gre_radial_golden_full_1071ms.seq'):
>>>>>>> 955e10bf0de467faa4ef33d7ec7fb494ae53de93
    # ======
    # SETUP
    # ======
     # FOV for SIEMENS prisma = 125 (3D), 250 (2D)
    fov = 250e-3 # initial val =260e-3
    Nx = 64  # Define FOV and resolution
    slice_thickness = 3e-3 #initial val = 3e-3  # Slice thickness
    n_slices = 1
    adc_oversampling = 2
    phi = np.pi/4

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
        adc_samples_limit=8192
    )

    seq = pp.Sequence(system)  # Create a new sequence object

    # ======
    # CREATE EVENTS
    # ======

    # Create fat saturation pulse

    B0 = 3
    sat_ppm = -3.35
    sat_freq=sat_ppm*1e-6*B0*system.gamma

    rf_fs = pp.make_gauss_pulse(
        flip_angle=110*(np.pi/180),
        system=system,
<<<<<<< HEAD
        duration=3e-3,
        dwell=10e-6,
        bandwidth=abs(sat_freq),
        
        apodization=0.5,
        time_bw_product=4,
=======
        duration=8e-3,
        dwell=10e-6,
        bandwidth=abs(sat_freq),
        freq_offset=sat_freq,
>>>>>>> 955e10bf0de467faa4ef33d7ec7fb494ae53de93
        use='saturation',
    )

    rf_fs.phase_ppm = -2*np.pi*rf_fs.freq_ppm*rf_fs.center_pos

    gz_fs = pp.make_trapezoid(
        channel='z',
        system=system,
        delay=pp.calc_duration(rf_fs),
        area=1/1e-4,
    )

    # Create 90-degree slice selection pulse and gradient
    rf, gz, gz_reph = pp.make_sinc_pulse(
        apodization=0.5,
        duration=3e-3,
        flip_angle=90 * np.pi / 180,
        slice_thickness=slice_thickness,
        system=system,
        time_bw_product=4,
        use='excitation',
        return_gz=True,

    )
   # Create raw single-shot archimedian spiral
   # define k-space parameters

    deltak = 1 / fov
    k_radius = np.round(Nx/2) #Spiral radius
    k_samples = np.round(2*np.pi*k_radius)*adc_oversampling #no. of samples on outest circle of spiral
    tos_calculation = 10 #time oversampling
    grad_oversampling = True
    c_max = k_radius*k_samples*tos_calculation
    r = np.array([])
    a = np.array([])
    ka = np.array([])

    for c in range(0, c_max):
        r_val = deltak*c/k_samples/tos_calculation
        r = np.append(r, r_val)

        a_val = (c % (k_samples*tos_calculation))*2*np.pi/k_samples/tos_calculation
        a = np.append(a, a_val)

        ka_val = r*np.exp(1j*a)
        ka = np.append(ka, ka_val)

    ka = [ka.real, ka.imag] #kx and ky matrix


    # Calc gradients and slew rates

    dt = system.grad_raster_time/tos_calculation
    ga, sa = pp.traj_to_grad(
        k=ka,
        raster_time=dt,
        first_grad_stephalf_raster=True,
        conservative_slew_est=True

    )

    # Limit analysis

    safety_margin = 0.99

    dt_gabs = abs(ga[0,:] + 1j*ga[1,:])/(system.max_grad*safety_margin)*dt #absolutely have no clue what this is and whether it will work
    dt_sabs = np.sqrt(abs(sa[0,:]+1j*sa[1,:]))/(system.max_slew*safety_margin)*dt

    dt_opt=np.max([dt_gabs, dt_sabs])

    #Apply lower limit to avoid losing trajectory detail

    dt_min = 4*system.grad_raster_time/k_samples/tos_calculation #minimum of 4 points per revolution
    dt_opt0 = dt_opt
    if dt_opt < dt_min:
        dt_opt = dt_min

    fig, ax = plt.figure()
    fig.subtitle('combined time stepping')
    ax.addsubplot(dt_opt0, dt_opt)
    fig.show()

    t_smooth = np.array([0,np.cumsum(dt_opt,2)] )

    dt_grad = system.grad_raster_time/(1+ grad_oversampling)

    if grad_oversampling:
        safety_factor_1st_timestep = 0.7
        t_end = t_smooth[-1] - (safety_factor_1st_timestep)*dt_grad
        #idk what this 0:np.floor() means?!! ARGH
        t_grad_x = [0]
        for x in range(0, (safety_factor_1st_timestep + (np.floor(t_end/dt_grad)))*dt_grad):
            t_grad = np.append(t_grad_x, x)

    else:
        t_end = t_smooth[-1] - 0.5*dt_grad
        t_grad_y = [0]
        for y in range(0, (safety_factor_1st_timestep + (np.floor(t_end/dt_grad)))*dt_grad):
            t_grad = np.append(t_grad_y,y)


    kopt = np.interp(np.conjugate((t_smooth, np.conjugate(ka), t_grad)))

    #analysis (ignored lines 86 to 96 for now of orig .m file)
    print("original duration" + str(round(1e6*dt*max(ka.shape))))
    print("duration smooth" + str(round(1e6*dt_grad*max(kopt.shape))))


    gos, sos = pp.traj_to_grad(
        k=kopt,
        raster_time=dt_grad,
        first_grad_stephalf_raster=grad_oversampling
        )
    spiral_grad_shape=gos
   #Ignored gradient and slew rate with abs constraint for now

    #Define gradients and adc events
    adc_time = dt_grad*spiral_grad_shape.shape[1]
    adc_samples_desired = k_radius*k_samples
    adc_dwell = round(adc_time/adc_samples_desired/system.adc_raster_time)*system.adc_raster_time
    adc_segments, adc_samples_persegment = pp.calc_adc_segments(
        num_samples=adc_samples_desired,
        dwell=adc_dwell,
        delay=round((pp.calc_duration(gz_reph)-adc_dwell/2)/system.rf_raster_time)*system.rf_raster_time
        )
    if grad_oversampling:
        spiral_grad_shape_initial = []
        for i in range
        #Ask a matlab bestie what is going on!!!!!



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
        with open('test_report_epi_fullsample.txt', 'w') as file:
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
<<<<<<< HEAD
        sar_values_table.to_csv('SAR_EPI_fullsample.csv', index=False)
=======
        sar_values_table.to_csv('SAR.csv', index=False)
>>>>>>> 955e10bf0de467faa4ef33d7ec7fb494ae53de93

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
<<<<<<< HEAD
    # ======
    # K-SPACE VISUALIZATION
    # ======
    if k_space:

        k_adc, k_all, t_excite, t_refocus, t_adc = seq.calculate_kspace()
        plt.figure()
        plt.title('Estimated K-Space Trajectory')
        plt.plot(k_adc[0], k_adc[1])
        plt.show()



=======
>>>>>>> 955e10bf0de467faa4ef33d7ec7fb494ae53de93



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
<<<<<<< HEAD
    main(plot=True, write_seq=True, pns_check=True, test_report=True, sar=False, acoustic_check=True, k_space=True)
=======
    main(plot=True, write_seq=True, sar=True)

>>>>>>> 955e10bf0de467faa4ef33d7ec7fb494ae53de93
