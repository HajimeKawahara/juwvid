module cohenclass
import jnufft

function tfrwv(x,y=NaN,t=NaN,f=NaN,N=NaN,silent=0,use_nufft=true)
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
        tau=round(N/2); 
        if ti<=xrow-tau && ti>=tau+1
            tfr[tau+1,icol] = 0.5*(x[ti+tau]*conj(y[ti-tau]) + x[ti-tau]*conj(y[ti+tau]))
        end
    end

    ###Choose FFT or DFT
    if isnan(f)[1]
        if silent==0 println("Use fft.") end
        for i=1:N
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    elseif use_nufft
        if silent==0 println("Use nufft.") end
        for i=1:N
            tfr[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-32)
        end
        return tfr
    else
        if silent==0 println("Use Direct DFT.") end
        Nf=size(f)[1]
        m=collect(1:N)
        tfrnew=zeros(Complex64,Nf,N)        
        for i=1:N
            for j=1:Nf                
                tfrnew[j,i]=sum(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            end            
        end
        return tfrnew
    end

    #--check--
    #for i=1:N
    #    tfr[:,i]=fft(tfr[:,i])
    #end
    #for i=1:N
    #    for j=1:Nf                
    #        println(j,"-",i," fft:",tfr[j,i]," dft:",tfrnew[j,i])
    #    end
    #end

end

function tfrpwv(x,y=NaN,t=NaN,f=NaN,N=NaN,h=NaN,silent=0)
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

    if isnan(f)[1]
        for i=1:N
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    else
        Nf=size(f)[1]
        m=collect(1:N)
        tfrnew=zeros(Complex64,Nf,N) 
        for i=1:N
            for j=1:Nf                
                tfrnew[j,i]=sum(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            end            
        end
        return tfrnew
    end

    return tfr
end

end

#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=cohenclass.tfrwv(z)
#nsample=16
#tfr=cohenclass.tfrwv(z,NaN,NaN,[1,2],NaN,0)
##tfr=cohenclass.tfrpwv(z)
#println(tfr)

