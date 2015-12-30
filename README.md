# juwvid

[![Licence](http://img.shields.io/badge/license-GPLv2-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

Julia codes for Cohen's class distribution. Currently, it's under development and functions available are very limited. I imported the Wigner-Ville distribution and the pseudo Wigner-Ville distribution from MATLAB GPL programs, tftb-0.2. Regarding tftb, visit http://tftb.nongnu.org/ .

## Requirement

- Julia v0.4

### Julia Packages 

- DSP
- Interpolations (in future)
- IJulia, Winston, Color (just for showing a sample)

```
Pkg.add("DSP")
```

## Available 

- Wigner-Ville distribution
- Cross Wigner-Ville distribution
- Pseudo Wigner-Ville distribution

## Planned 

- polynomial WV (possibly)

## Gallery

- Wigner-Ville

<img src="./figure/wv.png" Titie="explanation">

- Pseudo Wigner Ville

<img src="./figure/pwv.png" Titie="explanation">