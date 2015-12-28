module cohenclass
    #   These function were translated to this Julia code from MATLAB GPL programs, tftb-0.2.
    #   Licence is GPL v0.2 see License for the detail
    #   This program is free software; you can redistribute it and/or modify
    #   it under the terms of the GNU General Public License as published by
    #   the Free Software Foundation; either version 2 of the License, or
    #   (at your option) any later version.
    #    
    #   Original tftb MATLAB code was written by F. Auger, May-August 1994, July 1995.
    #	Regarding tftb, vist http://tftb.nongnu.org/ . 
    #
    #   Copyright (c) Hajime Kawahara (2015)

function tfrwv(x,t=NaN,N=NaN,silent=0)
    xcol,xrow = size(x) #xcol for cross-rwv
    if isnan(t) t=(1:xrow)' end
    if isnan(N) N=xrow end
    #information
#    println("Sizes of x and t")
#    println(size(x))
#    println(size(t))
    #error handling
    if N<0; println("N must be greater than zero"); exit(); end
    if xcol==0 || xcol>2; println("X must have one or two columns"); exit() end
    if xcol==1 && silent ==0 ; println("Single Wigner Ville"); end
    if xcol==2 && silent ==0 ; println("Cross Wigner Ville"); end
    if nextpow2(N)!=N && silent ==0 ; println("For a faster computation, N should be a power of two\n"); end

    tfr=zeros(Complex64,N,N) # plane by default

    for icol=1:N
        ti=t[icol]
        taumax=minimum([ti-1,xrow-ti,round(N/2)-1])
        tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = x[1,ti+tau].*conj(x[xcol,ti-tau]) #                
        tau=round(N/2); 
        if ti<=xrow-tau && ti>=tau+1
            tfr[tau+1,icol] = 0.5*(x[1,ti+tau]*conj(x[xcol,ti-tau]) + x[1,ti-tau]*conj(x[xcol,ti+tau]))
        end
    end
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
    if xcol==1
        tfr=real(tfr)
    end

    return tfr
end

function tfrpwv(x,t=NaN,N=NaN,h=NaN,silent=0)
    xcol,xrow = size(x) #xcol for cross-rwv
    if isnan(t) t=(1:xrow)' end
    if isnan(N) N=xrow end
    if isnan(h) 
        hlength=floor(N/4)
        hlength=hlength+1-rem(hlength,2)
        h=0.54 - 0.46*cos(2.0*pi*(1:hlength)'/(hlength+1)) #Hamming
    end
    
    hcol,hrow=size(h)
    Lh=round(Int,(hrow-1)/2)
    h=h/h[Lh+1]

    if hcol!=1 || rem(hrow,2)==0
        println("H must be a smoothing window with odd length")
        exit()
    end

    #information
    println("Sizes of x and t")
    println(size(x))
    println(size(t))

    #error handling
    if N<0; println("N must be greater than zero"); exit(); end
    if xcol==0 || xcol>2; println("X must have one or two columns"); exit() end
    if xcol==1 && silent ==0 ; println("pseudo Wigner Ville"); end
    if xcol==2 && silent ==0 ; println("Cross pseudo Wigner Ville"); end
    if nextpow2(N)!=N && silent ==0 ; println("For a faster computation, N should be a power of two\n"); end

    tfr=zeros(Complex64,N,N) # plane by default

    for icol=1:N
        ti=t[icol]
        taumax=minimum([ti-1,xrow-ti,round(N/2)-1,Lh])
        tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = h[Lh+1+tau]'.*x[1,ti+tau].*conj(x[xcol,ti-tau])
#        tfr[indices,icol] = x[1,ti+tau].*conj(x[xcol,ti-tau])

        tau=round(N/2); 
        if ti<=xrow-tau && ti>=tau+1 && tau<=Lh
            tfr[tau+1,icol] = 0.5*(h[Lh+1+tau]*x[1,ti+tau]*conj(x[xcol,ti-tau]) + h[Lh+1-tau]*x[1,ti-tau]*conj(x[xcol,ti+tau]))
        end
    end
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
    if xcol==1
        tfr=real(tfr)
    end

    return tfr
end

end

#import DSP
#y=[1.,2.,3.,5.,1.,2.,3.,5.,1.,2.,3.,5.,1.,2.,3.,5.]
#ya=DSP.Util.hilbert(y) # transpose is necessary 
#y=conj(ya')
#tfr=tfr4powv(y)
#tfr=tfrpwv(y)
#tfr=tfrwv(y)
#println(real(tfr))

