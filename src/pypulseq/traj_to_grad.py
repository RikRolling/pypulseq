from typing import Tuple, Union

import numpy as np

from pypulseq.opts import Opts


def traj_to_grad(k: np.ndarray, raster_time: Union[float, None] = None ,first_grad_stephalf_raster: bool = True, conservative_slew_estimate: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert k-space trajectory `k` into gradient waveform in compliance with `raster_time` gradient raster time.

    Parameters
    ----------
    k : numpy.ndarray
        K-space trajectory to be converted into gradient waveform.
    raster_time : float, default=Opts().grad_raster_time
        Gradient raster time.

    #RIK ADDITION

    first: numpy.ndarray

    first_grad_stephalf_raster: bool

    Conservative slew estimate: bool


    Returns
    -------
    g : numpy.ndarray
        Gradient waveform.
    sr : numpy.ndarray
        Slew rate.
    """
    if raster_time is None:
        raster_time = Opts.default.grad_raster_time

    # Compute finite difference for gradients in Hz/m
    g = (k[..., 1:] - k[..., :-1]) / raster_time
    # Compute the slew rate
    sr0 = (g[..., 1:] - g[..., :-1]) / raster_time
    if first_grad_stephalf_raster:
        sr0[...,0] = sr0[...,0]*2

    sr = np.zeros(sr0.shape[:-1] + (sr0.shape[-1] + 1,))
    sr[..., 0] = sr0[..., 0]

    if conservative_slew_estimate:
        if first_grad_stephalf_raster:
            sr[...,1] = sr0[...,1]
            sr[...,2:] = np.amax(sr0[...,1:-2], sr0[...,1:])
        else:
            sr[...,1:] = np.amax(sr0[:,1:-2], sr0[:,1:])
    else:
        if first_grad_stephalf_raster:
            sr[...,1]=sr0[...,1]
            sr[..., 2:] = 0.5 * (sr0[..., 1:-2] + sr0[..., 2:])
        else:

            sr[..., 1:] = 0.5 * (sr0[..., :-2] + sr0[..., 1:])
           # sr[..., -1] = sr0[..., -1]
















    #sr0[:,1] = sr0[:,1]*2


    # Gradient is now sampled between k-space points whilst the slew rate is between gradient points
    #sr = np.zeros(sr0.shape[:-1] + (sr0.shape[-1] + 1,))
    #sr[..., 0] = sr0[..., 0]
    #sr[..., 1:-1] = 0.5 * (sr0[..., :-1] + sr0[..., 1:])
    #sr[..., -1] = sr0[..., -1]

    return g, sr
