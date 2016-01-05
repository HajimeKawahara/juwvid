# juwvid

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

Julia codes for Cohen's class distribution. Currently, it's under development and functions available are very limited. I imported the Wigner-Ville distribution, the pseudo Wigner-Ville distribution, and the short-time Fourier transfrom from MATLAB GPL programs, tftb-0.2 and modified them. Regarding tftb, visit http://tftb.nongnu.org/ .

## Requirement

- Julia v0.4

### Julia Packages 

- DSP
- Interpolations (in future)
- IJulia, Winston, Color (just for showing a sample)

To install, use Pkg.add in the julia console, such as

```
Pkg.add("DSP")
```
## Julia Bindings of non-uniform FFT (NUFFT)

The Julia bindings of NUFFT are available (currently only for nufft1d2). The original fortran code of NUFFT is from [CMCL/NUFFT](http://www.cims.nyu.edu/cmcl/nufft/nufft.html) version 1.3.3. The code in the nufft_src directory is BSD licensed and the Julia bindings are GPLv2.

Use Makefile to generate a shared library file, libnufft.so. 

```
make
```



## Available 

- Wigner-Ville distribution (w/ NUFFT)
- Cross Wigner-Ville distribution  (w/ NUFFT)
- Pseudo Wigner-Ville distribution (w/ NUFFT)
- Short-time Fourier Transform
- polynomial WV (under development)




