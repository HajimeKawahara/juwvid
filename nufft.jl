#### DO NOT WORK YET !!!!!! ####

module nufft

function compute_grid_params(M, eps)
    # Choose Msp & tau from eps following Dutt & Rokhlin (1993)
    if eps <= 1E-33 or eps >= 1E-1
        println("eps =",eps,"must satisfy 1e-33 < eps < 1e-1.")
    end

    if eps > 1E-11 
        ratio = 2 
    else 
        ration = 3
    end
    
    Msp = round(Int64, -log(eps)/(pi*(ratio - 1)/(ratio - 0.5)) + 0.5)
    Mr = maximum(ratio*M, 2*Msp)
    lambda = Msp/(ratio*(ratio - 0.5))
    tau = pi*lambda/M^2
    
    return Msp, Mr, tau

end

function nufft1d(x, c, M, df=1.0, eps=1E-15, iflag=1)
    #"""Fast Non-Uniform Fourier Transform with Julia"""
    Msp, Mr, tau = compute_grid_params(M, eps)
    N = length(x)

    # Construct the convolved grid
    ftau = zeros(Complex64, Mr)
    Mr = size(ftau)[1]
    hx = 2*pi/Mr
    mm = round(Int64,rem(2*Msp+(-Msp:Msp-1),2*Msp)+1)

    for i in 1:N
        xi = (x[i] * df) % (2 * pi)
        m = 1 + round(Int, xi // hx)
        spread = exp(-0.25 * (xi - hx * (m + mm))^2 / tau)
        ftau[(m + mm) % Mr] += c[i] * spread

    # Compute the FFT on the convolved grid
    if iflag < 0
        Ftau = (1 / Mr) * fft(ftau)
    else
        Ftau = ifft(ftau)

    # irekae    
###    Ftau = np.concatenate([Ftau[-(M//2):], Ftau[:M//2 + M % 2]])

    # Deconvolve the grid using convolution theorem
    k = nufftfreqs(M)
    return (1/N)*sqrt(pi / tau) * exp(tau * k^2) * Ftau

end

end


########################################################
# Note
#
# nufft.jl was imported from nufftpy.
#
#Information on the licence of nufftpy: 
#------------------------------------------------------
#Copyright (c) 2014, Jake Vanderplas
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:

#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.

#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.