#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 08:23:53 2025
K-Space Tester

Converting image to k-space to explore undersampling
and image recon techniques
@author: ritikakhot
"""
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fftshift, ifftshift, fftn, ifftn


# Read image file
image = mpimg.imread('7T-MRI-Deep-Resolve-Brain.jpg')

plt.figure('Initial Image')
plt.imshow(image)
plt.title('Initial Image')
plt.show()


"""
Helpers for transforming data from k-space to image space and vice-versa.
from ismrmrd-pthon-tools
"""


def transform_kspace_to_image(k, dim=None, img_shape=None):
    """ Computes the Fourier transform from k-space to image space
    along a given or all dimensions

    :param k: k-space data
    :param dim: vector of dimensions to transform
    :param img_shape: desired shape of output image
    :returns: data in image space (along transformed dimensions)
    """
    if not dim:
        dim = range(k.ndim)

    img = fftshift(ifftn(ifftshift(k, axes=dim), s=img_shape, axes=dim), axes=dim)
    img *= np.sqrt(np.prod(np.take(img.shape, dim)))
    return img


def transform_image_to_kspace(img, dim=None, k_shape=None):
    """ Computes the Fourier transform from image space to k-space space
    along a given or all dimensions

    :param img: image space data
    :param dim: vector of dimensions to transform
    :param k_shape: desired shape of output k-space data
    :returns: data in k-space (along transformed dimensions)
    """
    if not dim:
        dim = range(img.ndim)

    k = fftshift(fftn(ifftshift(img, axes=dim), s=k_shape, axes=dim), axes=dim)
    k /= np.sqrt(np.prod(np.take(img.shape, dim)))
    return k


#Convert img to k-space

k_space = transform_image_to_kspace(image)
kx = k_space[0].real
ky = k_space[1].real
plt.figure('Initial K-Space')
plt.plot(kx, ky)
plt.show()






