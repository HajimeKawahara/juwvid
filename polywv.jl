module polywv
##### DO NOT WORKING YET ####
# Chandra Sekhar and Sreenivas 02
import Interpolations

function tfrpowv(x,t=NaN,N=NaN,silent=0)
    println("DON'T USE !!! DO NOT WORKING YET.")
    #f :frequency array
    # the six-order polynomial q=6
    xcol,xrow = size(x) #xcol for cross-rwv
    if isnan(t) t=(1:xrow) end
    if isnan(N) N=xrow end
    #information
#    println("Sizes of x and t")
#    println(size(x))
#    println(size(t))

    #error handling
    if N<0; println("N must be greater than zero"); exit(); end
    if xcol==0 || xcol>2; println("X must have one or two columns"); exit() end
    if xcol==1 && silent==0; println("Single q=6 polynomial Wigner Ville"); end
    if xcol==2 && silent==0; println("Cross  q=6 polynomial Wigner Ville"); end

    d1=0.675
    d2=0.85
    #interpolation
    za = Interpolations.interpolate((t,),conj(x[1,:]'[:,1]), Interpolations.Gridded(Interpolations.Linear()));
    zb = Interpolations.interpolate((t,),conj(x[xcol,:]'[:,1]), Interpolations.Gridded(Interpolations.Linear()));

    tau = 1/0.675
#    tau = 0.5/0.675
    M=N
    tfr=zeros(Complex64,M,M) # plane by default

    for icol=1:N
        ti=t[icol]
        println("---",icol)
        for mrow=-M+1:M-1                        
            if ti+d2*mrow*tau <= t[end] && ti-d2*mrow*tau >= 1                
                println(mrow)
                tfr[mrow,icol] = za[ti+d1*mrow*tau].*zb[ti+d1*mrow*tau].*conj(zb[ti-d1*mrow*tau]).*conj(zb[ti-d1*mrow*tau]).*conj(zb[ti+d2*mrow*tau]).*zb[ti-d2*mrow*tau] #                
            end
        end
    end

    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
#    println(tfr)
#    exit()
    if xcol==1
        tfr=real(tfr)
    end

    return tfr
end
end

import DSP
y=linspace(0.0,8.0,8)
ya=DSP.Util.hilbert(y) # transpose is necessary 
y=conj(ya')
tfr=polywv.tfrpowv(y)
