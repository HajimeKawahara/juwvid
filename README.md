# juwvid

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/polydata.png" Titie="explanation" Width=200px>
<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/polywv1.png" Titie="explanation" Width=200px>
<img src="https://github.com/HajimeKawahara/juwvid/blob/master/figure/atfr.png" Titie="explanation" Width=200px>

Julia codes for time-frequency analysis. I started to import the Wigner distribution, the pseudo Wigner distribution, and the short-time Fourier transform from MATLAB GPL programs, tftb-0.2 and have modified them. Regarding tftb, visit http://tftb.nongnu.org/. The modification includes the non-uniform FFT, the adaptive algorithm by Stankovic, Dakovic, Thayaparan 2013, the S-method, the L-Wigner distribution, and the polynomial Wigner-Ville distribution.

#### References 
- Cohen, L. 1995, Time-Frequency Analysis (PTR-PH)
- Stankovic, L., Dakovi ́c, M., & Thayaparan, T. 2013, Time-frequency signal analysis with applications (Artech House)
- Boashash, B. 2015, Time-Frequency Signal Analysis and Processing, 2nd Edition A Comprehensive Reference (Elsevier)

## Requirement

- Julia v0.4

#### Julia Packages 

- DSP
- Jupyter, Winston, Color, PyPlot (these are just for tutorials)

To install, use Pkg.add in the julia console, such as

```
Pkg.add("DSP")
```

#### Julia Bindings of non-uniform FFT (NUFFT)

Juwvid has the Julia bindings of parts of [CMCL/NUFFT](http://www.cims.nyu.edu/cmcl/nufft/nufft.html) version 1.3.3 (but currently only for nufft1d2). The original Fortran code of NUFFT is from [CMCL/NUFFT](http://www.cims.nyu.edu/cmcl/nufft/nufft.html). The source files in the nufft_src directory are BSD licensed, and the Julia bindings are GPLv2.

Edit Makefile and type

```
make
```

to generate a shared library file, libnufft.so. Change the path of ccall in jnufft.jl. 

## Available 

- Wigner-Ville distribution w/ NUFFT
- Cross Wigner-Ville distribution w/ NUFFT
- Pseudo Wigner-Ville distribution w/ NUFFT and Adaptive Algorithm
- Short-time Fourier Transform w/ NUFFT
- S-method w/ NUFFT
- L-Wigner distribution w/ NUFFT
- polynomial Wigner-Ville distribution w/ NUFFT

#### Examples
See ipython notebooks

