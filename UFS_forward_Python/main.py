import numpy as np
import matplotlib.pyplot as plt

from forward import patterns, params2fairSIM, psf, forward

N = 1024
sz = (N, N)  # Size of superresolved image
num_phases = 3
num_orientations = 3
num_patterns = num_phases * num_orientations
num_t = num_patterns * 2  # Number of time points = number of acquired SIM images (multiple of number of patterns)

# Generation of theoretical patterns based on physical parameters
phase = np.linspace(0, 2*np.pi, num_phases, endpoint=False)
orientations = np.array([0, np.pi/3, 2*np.pi/3])
Na = 1.49  # Numerical aperture
nl = 1.515  # Refractive index of the objective medium (glass/oil)
ns = 1.333  # Refractive index of the sample medium (water)
lamb = 605  # Emission wavelength (in nm)
res = 32  # Resolution of superresolved image (in nm)
# These are the parameters estimated by FairSIM
phase_offset, shift = params2fairSIM(orientations, phase, Na, ns, nl, lamb, res, sz)
patterns = patterns(sz, num_phases, phase_offset, shift, a=1)

plt.figure(1)
plt.imshow(patterns[:, :, 0])
plt.colorbar()

# Simulated PSF based on physical parameters
psf, otf = psf(sz, lamb, res, Na)

# Example of application of forward operator
vol = np.random.rand(N, N, num_t)
y = forward(vol, otf, patterns)

plt.figure(2)
plt.imshow(y[:, :, 0])
plt.colorbar()
plt.show()

