module cohenclass

function tfrwv(x,y=NaN,t=NaN,N=NaN,silent=0)
    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end
    if isnan(N) N=xrow end
    if isnan(y)[1]
        if silent ==0  println("Single Wigner Ville") end
        y=x
        sw=0
    else 
        if silent==0 println("Cross Wigner Ville") end
        sw=1
    end
    tfr=zeros(Complex64,N,N) # plane by default
    for icol=1:N
        ti=t[icol]
        taumax=minimum([ti-1,xrow-ti,round(N/2)-1])
        tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = x[ti+tau].*conj(y[ti-tau]) #                
#        for i in 1:length(indices)
#            ii=indices[i]
#            println(ii,"   ",x[ti+tau[i]].*conj(y[ti-tau[i]]))
#        end
#        println(fft(tfr[:,icol]))
        tau=round(N/2); 
        if ti<=xrow-tau && ti>=tau+1
            tfr[tau+1,icol] = 0.5*(x[ti+tau]*conj(y[ti-tau]) + x[ti-tau]*conj(y[ti+tau]))
        end
    end
    ############
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
    if sw==0
        tfr=real(tfr)
    end

    return tfr
end

function tfrpwv(x,y=NaN,t=NaN,N=NaN,h=NaN,silent=0)
    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end
    if isnan(N) N=xrow end
    if isnan(y)[1] 
        if silent ==0 println("Single pseudo Wigner Ville") end
        y=x
        sw=0
    else
        if silent ==0 println("Cross pseudo Wigner Ville") end
        sw=1
    end
    if isnan(h)[1] 
        hlength=floor(N/4)
        hlength=hlength+1-rem(hlength,2)
        h=0.54 - 0.46*cos(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
    end    
    hrow = size(h)[1] 
    Lh=round(Int,(hrow-1)/2)
    h=h/h[Lh+1]

    tfr=zeros(Complex64,N,N) # plane by default

    for icol=1:N
        ti=t[icol]
        taumax=minimum([ti-1,xrow-ti,round(N/2)-1,Lh])
        tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = h[Lh+1+tau].*x[ti+tau].*conj(y[ti-tau])
        tau=round(N/2); 
        if ti<=xrow-tau && ti>=tau+1 && tau<=Lh
            tfr[tau+1,icol] = 0.5*(h[Lh+1+tau]*x[ti+tau]*conj(y[ti-tau]) + h[Lh+1-tau]*x[ti-tau]*conj(y[ti+tau]))
        end
    end
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
    if sw==0
        tfr=real(tfr)
    end

    return tfr
end

end

#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=cohenclass.tfrwv(z)
##tfr=cohenclass.tfrpwv(z)
#println(real(tfr))

