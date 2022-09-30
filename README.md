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
- The folder **Utils** contains auxiliary functions.
- The folder **UFS_forward_Python** contains Python code implementing the time-varying SIM forward model. The function that generate the illumination patterns takes parameters estimated by FairSIM as input.

## Additional comments / things to look into

- The script **Utils/pattern_parameters_FairSIM.java** allows to estimate the parameters of the illumination patterns using FairSIM (ImageJ plugin) and write them in a .txt file, which is then parsed by **PatternsFromFairSimTxt.m** within the **RealDataUFS.m** script. Note that playing with the OTF dampening parameter **otfCorr** can significantly affect the parameter estimation. Currently, the estimated parameters are simply averaged over all the 9-frame cycles. With the BIOP data I used, the variability in the estimated parameters is quite low among cycles, but the quality of the estimation is still considered as "weak" by FairSIM.

- FairSIM does not estimate the amplitude parameter **a** of the illumination patterns (1 + a * cos(...)). I tried using Alejandro's FlexSIM method to estimate the pattern parameters (including the amplitude), I didn't have time to investigate much but when I ran it as is there is a huge variability in the estimation across 9-frame cycles, so it would require more work. Note that it estimated the amplitude parameter around a=0.01, which may explain why the FairSIM estimation is weak. The is consistent with the fact that the patterns are very hard to see "manually" (even in Fourier domain), so the acquisition may have been a bit crappy (Arne had kind of acknowledged that, but I first tried with what we had).

- Concerning the OTF, I generate an estimate based on acquisition parameters using Manu's code, I have no clue how accurate it is. FairSIM and Alejandro both use a different estimation of the OTF, with a dampening parameter but without the refractive indices (I didn't check how they compare).
