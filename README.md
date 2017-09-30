# juwvid

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![arXiv](http://img.shields.io/badge/arXiv-1603.02898-green.svg?style=flat)](http://arxiv.org/abs/1603.02898)
[![ascl](http://img.shields.io/badge/ascl-1702.003-red.svg?style=flat)](http://ascl.net/1702.003)

<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/nufft.png" Titie="explanation" Width=400px>
<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/smnufft.png" Titie="explanation" Width=200px>
<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/ligopwv.png" Titie="explanation" Width=400px>

Julia codes for time-frequency analysis. I started to import the Wigner distribution, the pseudo Wigner distribution, and the short-time Fourier transform from MATLAB GPL programs, tftb-0.2 and have modified them. Regarding tftb, visit http://tftb.nongnu.org/. The modification includes the zero-padding FFT, the non-uniform FFT, the adaptive algorithm by Stankovic, Dakovic, Thayaparan 2013, the S-method, the L-Wigner distribution, and the polynomial Wigner-Ville distribution. Juwvid was originally developed for my paper, [Kawahara (2016)](http://arxiv.org/abs/1603.02898).

## Available 

- Wigner-Ville distribution w/ NUFFT
- Cross Wigner-Ville distribution w/ NUFFT
- Pseudo Wigner-Ville distribution w/ zero-padding FFT and NUFFT and Adaptive Algorithm
- Short-time Fourier Transform w/ NUFFT
- S-method w/ NUFFT
- L-Wigner distribution w/ NUFFT
- polynomial Wigner-Ville distribution w/ NUFFT

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

## Tutorials

See ipython notebooks in the ipynb directory.

- Wigner Ville Distribution with FFT, zero-padding FFT, and NUFFT.ipynb
- Pseudo Wigner Ville Distribution for the nonlinear IF using FFT and NUFFT.ipynb
- Adaptive algorithm for window selection of the pseudo WV.ipynb
- S-Method.ipynb
- S-Method with NuFFT.ipynb
- L-Wigner Distribution and Polynomial Wigner Ville Distribution.ipynb
- Polynomial Wigner Ville distribution with NuFFT.ipynb
- Application to LIGO data.ipynb [PDF](https://github.com/HajimeKawahara/juwvid/blob/master/documents/Application\ to\ LIGO\ data.pdf)

#### References 
- Cohen, L. 1995, Time-Frequency Analysis (PTR-PH)
- Stankovic, L., Dakovi ́c, M., & Thayaparan, T. 2013, Time-frequency signal analysis with applications (Artech House)
- Boashash, B. 2015, Time-Frequency Signal Analysis and Processing, 2nd Edition A Comprehensive Reference (Elsevier)

