# juwvid

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

Julia codes for Cohen's class distribution. Currently, it's under development and functions available are very limited. I imported the Wigner-Ville distribution, the pseudo Wigner-Ville distribution, and the short-time Fourier transform from MATLAB GPL programs, tftb-0.2 and modified them. Regarding tftb, visit http://tftb.nongnu.org/. The modification includes the non-uniform FFT and the adaptive algorithm by Stankovic, Dakovic, Thayaparan 2013. 

## Requirement

- Julia v0.4

#### Julia Packages 

- DSP
- Interpolations (in future)
- IJulia, Winston, Color (just for showing a sample)

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

- Wigner-Ville distribution (w/ NUFFT)
- Cross Wigner-Ville distribution  (w/ NUFFT)
- Pseudo Wigner-Ville distribution (w/ NUFFT and Adaptive Algorithm)
- Short-time Fourier Transform

#### Examples
See ipython notebooks

