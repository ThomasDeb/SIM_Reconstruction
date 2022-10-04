import numpy as np
import scipy as sp
from scipy import signal


def patterns(sz, num_phases, phase_offsets, shifts, a=1):
    """
    Generate patterns based on FairSIM data
    :param sz: (tuple) size of reconstructed image (twice the resolution of the data)
    :param num_phases: (int) number of phase shifts per orientation
    :param phase_offsets: (num_orientations) phase offsets for each orientation (relative to desired phase shift)
    :param shifts: (2 x num_orientations) orientation vectors
    :param a: (float) amplitude coefficient of patterns (not estimated by FairSIM)
    :return: patterns: (sz x num_patt) patterns
    """
    num_orientations = shifts.shape[1]
    num_patt = num_phases * num_orientations
    patt = np.zeros(sz + (num_patt,))
    phase_shifts = np.linspace(0, 2 * np.pi, num_phases, endpoint=False)

    k0 = np.linalg.norm(shifts, axis=0)
    orientations = - np.arctan2(shifts[1, :], shifts[0, :])  # Minus is due to a FairSIM convention?

    x, y = np.meshgrid(np.arange(0, sz[0]), np.arange(0, sz[1]))
    p = 0
    for i in range(num_orientations):
        k = 2 * np.pi * np.array([np.cos(orientations[i]), np.sin(orientations[i])]) * k0[i] / (2 * np.array(sz))
        for j in range(num_phases):
            patt[:, :, p] = 1 + a * np.cos(2 * (x * k[0] + y * k[1]) + phase_offsets[i] + phase_shifts[j])
            p += 1

    return patt


def params2fairSIM(orientations, phases, Na, ns, nl, lamb, res, sz):
    """
    Convert physical parameters to FairSIM parameters to generate patterns
    :param orientations: (num_orientations) orientation angles
    :param phases: (num_phases) list of phases for each orientation
    :param res: (float) resolution of a pixel of superresolved image (in nm)
    :param Na: (float) numerical aperture
    :param ns:
    :param nl:
    :param lamb: (float) emission wavelength (in nm)
    :param sz: (tuple) size of reconstructed image (twice the resolution of the data)
    :return:
    """
    k0 = Na * res * ns / (nl * lamb)  # Obtained by reverse engineering Manu's code (dimensionless)
    shifts = np.stack([np.cos(-orientations), np.sin(-orientations)], axis=0) * k0 * 2 * np.reshape(sz, (2, 1))
    # Minus due to FairSIM convention
    phase_shifts = np.linspace(0, 2*np.pi, len(phases), endpoint=False)  # Desired phase shifts
    phase_offsets = phases - phase_shifts  # Offsets with respect to phase
    return phase_offsets, shifts


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return rho, theta


def psf(sz, lamb, res, Na):
    fc = 2 * Na / lamb * res
    lx = np.linspace(-0.5, 0, int(sz[0] / 2), endpoint=False)
    rx = np.linspace(0, 0.5, int(sz[0] / 2))
    gridx = np.concatenate((lx, rx))  # To have a grid that includes 0
    ly = np.linspace(-0.5, 0, int(sz[1] / 2), endpoint=False)
    ry = np.linspace(0, 0.5, int(sz[1] / 2))
    gridy = np.concatenate((ly, ry))  # To have a grid that includes 0

    x, y = np.meshgrid(gridx, gridy)
    rho, theta = cart2pol(x, y)
    temp = 2 * np.arccos(np.minimum(rho / fc, 1))  # Argument is 1 in cutoff zone -> arccos is then 0
    otf = np.fft.fftshift(1 / np.pi * (temp - np.sin(temp)) * (rho < fc))
    psf = np.real(np.fft.fftshift(np.fft.ifft2(otf)))  # Ensure that psf is real
    otf = np.fft.fft2(np.fft.fftshift(psf))
    return psf, otf


def conv2D(im, otf):
    return np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(np.fft.fftshift(im)) * otf)))


def conv_volume(vol, otf):
    if vol.ndim == 2:
        out = conv2D(vol, otf)
    else:
        out = np.zeros_like(vol)
        for i in range(vol.shape[2]):
            out[:, :, i] = conv2D(vol[:, :, i], otf)
    return out


def forward(vol, otf, patterns):
    num_t = vol.shape[2]
    num_patt = patterns.shape[2]
    y = np.zeros((int(vol.shape[0] / 2), int(vol.shape[1] / 2), num_t))
    for i in range(int(num_t / num_patt)):
        rng = range(num_patt * i, num_patt * (i+1))
        y[:, :, rng] = conv_volume(vol[:, :, rng] * patterns, otf)[::2, ::2, :]
    return y
