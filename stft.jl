module stft
#   These functions were imported from MATLAB GPL programs,tftb-0.2
#   and were modified (simplified) by HK.
#   Licence is GPL v2 see License for the detail
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#    
#   Original tftb MATLAB code was written by F. Auger, May-August 1994, July 1995.
#	Regarding tftb, vist http://tftb.nongnu.org/ . 
#
#   Copyright (c) Hajime Kawahara (2015)

function tfrstft(x,t=NaN,N=NaN,h=NaN)
    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end
    if isnan(N) N=xrow end
    if isnan(h)[1] 
        hlength=floor(N/4)
        hlength=hlength+1-rem(hlength,2)
        h=0.54 - 0.46*cos(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
    end        
    h=h/norm(h)
    hrow=length(h)
    Lh=round(Int64,(hrow-1)/2)
 
    tfr=zeros(Complex64,N,N) # plane by default
    for icol=1:N
        ti=t[icol]
     # tau=collect(-minimum([round(N/2)-1,Lh,ti-1]):minimum([round(N/2)-1,Lh,xrow-ti]))
     # indices= rem(N+tau,N)+1; 
     # tfr(indices,icol)=x(ti+tau,1).*conj(h(Lh+1+tau));

        taumin=-minimum([round(N/2)-1,Lh,ti-1])
        taumax=minimum([round(N/2)-1,Lh,xrow-ti])
        tau=round(Int64,taumin:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = x[ti+tau].*conj(h[Lh+1+tau])
    end

#    tfr=fft(tfr)
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end

return tfr

end
end
#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=stft.tfrstft(z)
##tfr=cohenclass.tfrpwv(z)
#println(tfr)
