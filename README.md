<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/juwvid.png" Titie="explanation" Width=200px>

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![arXiv](http://img.shields.io/badge/arXiv-1603.02898-green.svg?style=flat)](http://arxiv.org/abs/1603.02898)
[![arXiv](http://img.shields.io/badge/arXiv-1810.00334-green.svg?style=flat)](http://arxiv.org/abs/1810.00334)
[![ascl](http://img.shields.io/badge/ascl-1702.003-red.svg?style=flat)](http://ascl.net/1702.003)

Julia codes for time-frequency analysis. I started to import the Wigner distribution, the pseudo Wigner distribution, and the short-time Fourier transform from MATLAB GPL programs, tftb-0.2 and modified them and added other techniques. Regarding tftb, visit http://tftb.nongnu.org/. The modification and new additions include the zero-padding FFT, the non-uniform FFT, the adaptive algorithm by Stankovic, Dakovic, Thayaparan 2013, the S-method, the L-Wigner distribution, the polynomial Wigner-Ville distribution, extraction of instantaneous frequency, mode tracking, time-frequency distributions of the Stokes parameters for the spectrogram, the Wigner-Ville distribution, the S-method. Juwvid was originally developed for the paper on a characterization method of Earth-like exoplanets, [Kawahara (2016)](http://arxiv.org/abs/1603.02898) (TFA in general). It was extended for the paper on the gravitational wave of core collapse super novae [Kawahara et al. (2018)](http://arxiv.org/abs/1810.00334) (polarization and IF tracking). If you want to refer juwvid, consider to cite those papers because the detailed description are ginen in them, depending on topics though.

## Current Stable Version

[Juwvid 0.6](https://github.com/HajimeKawahara/juwvid/releases/tag/v0.6)

## Available 

- Wigner-Ville distribution w/ NUFFT
- Cross Wigner-Ville distribution w/ NUFFT
- Pseudo Wigner-Ville distribution w/ zero-padding FFT and NUFFT and Adaptive Algorithm
- Short-time Fourier Transform w/ NUFFT
- S-method w/ NUFFT
- L-Wigner distribution w/ NUFFT
- polynomial Wigner-Ville distribution w/ NUFFT
- Stokes parameters for the Short-time Fourier Transform (spoectrogram-type Stokes)
- Stokes parameters for the Wigner-Ville distribution
- Stokes parameters for the pseudo Wigner-Ville distribution
- Stokes parameters for the S-method (S-type Stokes)

## Requirements and Install

- Julia v0.6 (VERSION<0.6 does not work now.)

#### Julia Packages 

- DSP
- PyPlot
- Distributions, IJulia (just for tutorials)

To install, use Pkg.add in the julia console, such as

```
Pkg.add("DSP")
```

Set JULIA_LOAD_PATH as 
```
setenv JULIA_LOAD_PATH /Users/kawahara/juwvid
```
for csh/tcsh or 
```
export JULIA_LOAD_PATH=/Users/kawahara/juwvid
```
for bash. 

#### Julia Bindings of non-uniform FFT (NUFFT)

Juwvid has the Julia bindings of parts of [CMCL/NUFFT](http://www.cims.nyu.edu/cmcl/nufft/nufft.html) version 1.3.3 (but currently only for nufft1d2). The original Fortran code of NUFFT is from [CMCL/NUFFT](http://www.cims.nyu.edu/cmcl/nufft/nufft.html). The source files in the nufft_src directory are BSD licensed, and the Julia bindings are GPLv2.

Edit Makefile and type

```
make
```

to generate a shared library file, libnufft.so. Change the path of ccall in jnufft.jl. 

#### Features

- Various TFDs
- Dense frequency sampling using NUFFT
- Thinning out the time grid for large dataset
- TFDs for polarization

## Tutorials

See ipython notebooks in the ipynb directory.

- Wigner Ville Distribution with FFT, zero-padding FFT, and NUFFT.ipynb
- Pseudo Wigner Ville Distribution for the nonlinear IF using FFT and NUFFT.ipynb
- Adaptive algorithm for window selection of the pseudo WV.ipynb
- S-Method.ipynb
- S-Method with NuFFT.ipynb
- L-Wigner Distribution and Polynomial Wigner Ville Distribution.ipynb
- Polynomial Wigner Ville distribution with NuFFT.ipynb
- Stokes Distribution.ipynb
- Application to LIGO data.ipynb 


#### References 
- Cohen, L. 1995, Time-Frequency Analysis (PTR-PH)
- Stankovic, L., Dakovi ́c, M., & Thayaparan, T. 2013, Time-frequency signal analysis with applications (Artech House)
- Boashash, B. 2015, Time-Frequency Signal Analysis and Processing, 2nd Edition A Comprehensive Reference (Elsevier)
- Kawahara, H. 2016, Frequency Modulation of Directly Imaged Exoplanets: Geometric Effect as a Probe of Planetary Obliquity, ApJ 822, 112 
- Kawahara, H. et al. 2018, A Linear and Quadratic Time-Frequency Analysis of Gravitational Waves from Core-Collapse Supernovae, accepted for publication in ApJ

