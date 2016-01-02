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

## Available 

- Wigner-Ville distribution
- Cross Wigner-Ville distribution
- Pseudo Wigner-Ville distribution
- Short-time Fourier Transform
- polynomial WV (under development)

## Sample

See https://gist.github.com/HajimeKawahara/398c99cdcd9863c8a8d9

