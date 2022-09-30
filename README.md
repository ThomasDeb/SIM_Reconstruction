# SIM_Reconstruction

## Description

Time-varying 2D structured illumination microscopy (SIM) reconstruction algorithm based on the inner-loop-free ADMM algorithm from [1]. As in [1], we use Hessian-Schatten regularization in the spatial dimensions. We use an additional total-variation (TV) regularizer for the temporal dimension, either with periodic boundary conditions which allows to adapt the inner-loop-free method of [1], or with zero boundary conditions using the direct proximity operator of 1D TV.

[1] <a href="https://ieeexplore.ieee.org/document/8579117" target="_blank">Computational Super-Sectioning for Single-Slice Structured Illumination Microscopy.</a>, <br />
IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019.  <br />
E. Soubies, and M. Unser.

[2] <a href="https://ieeexplore.ieee.org/document/6579659" target="_blank">A Direct Algorithm for 1-D Total Variation Denoising</a>, <br />
IEEE Signal Processing Letters, vol. 20 no. 11, pp. 1054-1057, 2013.  <br />
L. Condat.

## Requirements

The code requires the GlobalBioIm library v1.1.2 (or more recent releases) <br />
https://biomedical-imaging-group.github.io/GlobalBioIm/

## Repository content

The repository is organized as follows.

- The script **SimuUFS_MT.m** contains the main code for our method on simulated data.
- The script **SimuUFS_MT.m** allows to generate simulated time-varying ground-truth microtubule-like structures and to generate simulated acquisition parameters such as OTF and illumination patterns (used in SimuUFS_MT.m).
- The script **main_UFS_RealData.m** contains the main code for our method on real data.
- The script **RealDataUFS.m** to generate estimated acquisition parameters such as OTF and illumination patterns (used in main_UFS_RealData.m).
- Folder **Utils** contains auxilliary functions.
