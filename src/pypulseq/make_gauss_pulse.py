import math
from copy import copy
from types import SimpleNamespace
from typing import Tuple, Union
from warnings import warn

import numpy as np

from pypulseq.make_trapezoid import make_trapezoid
from pypulseq.opts import Opts
from pypulseq.supported_labels_rf_use import get_supported_rf_uses
from pypulseq.utils.tracing import trace, trace_enabled


def make_gauss_pulse(
    flip_angle: float,
    apodization: float = 0,
    bandwidth: float = 0,
    center_pos: float = 0.5,
    delay: float = 0,
    dwell: float = 0,
    duration: float = 4e-3,
    freq_offset: float = 0,
    max_grad: float = 0,
    max_slew: float = 0,
    phase_offset: float = 0,
    return_gz: bool = False,
    slice_thickness: float = 0,
    system: Union[Opts, None] = None,
    time_bw_product: float = 4,
    use: str = str(),
) -> Union[
    SimpleNamespace,
    Tuple[SimpleNamespace, SimpleNamespace, SimpleNamespace],
]:
    """
    Create a [optionally slice selective] Gauss pulse.
    RIK TO ADD FREQPPM AND PHASEPPM AS PER MATLAB FUNCTIONALITY!
    See also `pypulseq.Sequence.sequence.Sequence.add_block()`.

    Parameters
    ----------
    flip_angle : float
        Flip angle in radians.
    apodization : float, default=0
        Apodization.
    bandwidth : float, default=0
        Bandwidth in Hertz (Hz).
    center_pos : float, default=0.5
        Position of peak.
    delay : float, default=0
        Delay in seconds (s).
    dwell : float, default=0
    duration : float, default=4e-3
        Duration in seconds (s).
    freq_offset : float, default=0
        Frequency offset in Hertz (Hz).
    max_grad : float, default=0
        Maximum gradient strength of accompanying slice select trapezoidal event.
    max_slew : float, default=0
        Maximum slew rate of accompanying slice select trapezoidal event.
    phase_offset : float, default=0
        Phase offset in Hertz (Hz).
    return_gz : bool, default=False
        Boolean flag to indicate if the slice-selective gradient has to be returned.
    slice_thickness : float, default=0
        Slice thickness of accompanying slice select trapezoidal event. The slice thickness determines the area of the
        slice select event.
    system : Opts, default=Opts()
        System limits.
    time_bw_product : int, default=4
        Time-bandwidth product.
    use : str, default=str()
        Use of radio-frequency gauss pulse event. Must be one defined in pypulseq.supported_labels_rf_use.get_supported_rf_uses.

    Returns
    -------
    rf : SimpleNamespace
        Radio-frequency gauss pulse event.
    gz : SimpleNamespace, optional
        Accompanying slice select trapezoidal gradient event.
    gzr : SimpleNamespace, optional
        Accompanying slice select rephasing trapezoidal gradient event.

    Raises
    ------
    ValueError
        If invalid `use` is passed.
        If `return_gz=True` and `slice_thickness` was not passed.
    """
    if system is None:
        system = Opts.default

    if use != '' and use not in get_supported_rf_uses():
        raise ValueError(f'Invalid use parameter. Must be one of {get_supported_rf_uses()}. Passed: {use}')

    if dwell == 0:
        dwell = system.rf_raster_time

    if bandwidth == 0:
        bandwidth = time_bw_product / duration
    else:
        bandwidth = bandwidth
    alpha = apodization
    n_samples = round(duration / dwell)
    t = (np.arange(1, n_samples + 1) - 0.5) * dwell
    tt = t - (duration * center_pos)
    window = 1 - alpha + alpha * np.cos(2 * np.pi * tt / duration)
    signal = window * __gauss(bandwidth * tt)
    flip = np.sum(signal) * dwell * 2 * np.pi
    signal = signal * flip_angle / flip

    rf = SimpleNamespace()
    rf.type = 'rf'
    rf.signal = signal
    rf.t = t
    rf.shape_dur = n_samples * dwell
    rf.freq_offset = freq_offset
    rf.phase_offset = phase_offset
    rf.dead_time = system.rf_dead_time
    rf.ringdown_time = system.rf_ringdown_time
    rf.delay = delay
    if use != '':
        rf.use = use

    if rf.dead_time > rf.delay:
        warn(
            f'Specified RF delay {rf.delay * 1e6:.2f} us is less than the dead time {rf.dead_time * 1e6:.0f} us. Delay was increased to the dead time.',
            stacklevel=2,
        )
        rf.delay = rf.dead_time
    # RK added extra line
    rf.center=duration*center_pos

    if return_gz:
        if slice_thickness == 0:
            raise ValueError('Slice thickness must be provided')

        if max_grad > 0:
            system = copy(system)
            system.max_grad = max_grad

        if max_slew > 0:
            system = copy(system)
            system.max_slew = max_slew

        amplitude = bandwidth / slice_thickness
        area = amplitude * duration
        gz = make_trapezoid(channel='z', system=system, flat_time=duration, flat_area=area)
        gzr = make_trapezoid(
            channel='z',
            system=system,
            area=-area * (1 - center_pos) - 0.5 * (gz.area - area),
        )

        if rf.delay > gz.rise_time:
            gz.delay = math.ceil((rf.delay - gz.rise_time) / system.grad_raster_time) * system.grad_raster_time

        if rf.delay < (gz.rise_time + gz.delay):
            rf.delay = gz.rise_time + gz.delay

    # Following 2 lines of code are workarounds for numpy returning 3.14... for np.angle(-0.00...)
    negative_zero_indices = np.where(rf.signal == -0.0)
    rf.signal[negative_zero_indices] = 0

    if trace_enabled():
        rf.trace = trace()

    if return_gz:
        return rf, gz, gzr
    else:
        return rf


def __gauss(x: np.ndarray) -> np.ndarray:
    return np.exp(-np.pi * np.square(x))
