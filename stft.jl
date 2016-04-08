module stft
import jnufft

function tfrstft(x,t=NaN,N=NaN,f=NaN,itc=NaN,h=NaN,nwindow=4,silent=0,use_nufft=true)
    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end
    if isnan(N) N=xrow end

    if isnan(itc)[1]  
        Nt=N
        itc=collect(1:Nt)
    else
        Nt=length(itc)
    end

    if isnan(h)[1] 
        hlength=floor(N/nwindow)
#        hlength=floor(N/4)
        hlength=hlength+1-rem(hlength,2)
        h=0.54 - 0.46*cos(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
    end        
    h=h/norm(h)
    hrow=length(h)
    Lh=round(Int64,(hrow-1)/2) ##??

    tfr=zeros(Complex64,N,Nt) # plane by default
    for icol=1:Nt       
#        ti=t[icol]
        ti=t[itc[icol]]
        taumin=-minimum([floor(N/2)-1,Lh,ti-1])
        taumax=minimum([floor(N/2)-1,Lh,xrow-ti])
        tau=round(Int64,taumin:taumax)
        indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = x[ti+tau].*conj(h[Lh+1+tau])
    end

#    tfr=fft(tfr)
    ###Choose FFT or DFT
    if isnan(f)[1] 
        if silent==0 println("Use fft.") end
        for i=1:Nt
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    elseif use_nufft 
        if silent==0 println("Use nufft.") end
        Nf=size(f)[1]
        tfrnew=zeros(Complex64,Nf,N) 
        for i=1:Nt
            tfrnew[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-28)[1:Nf]
        end        
        return tfrnew
    else
        if silent==0 println("Use Direct DFT.") end
        Nf=size(f)[1]
        m=collect(1:Nt)
        tfrnew=zeros(Complex64,Nf,N)        
        for i=1:Nt
            for j=1:Nf                
                tfrnew[j,i]=sum(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            end            
        end
        return tfrnew
    end

end
end
#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=stft.tfrstft(z)
##tfr=cohenclass.tfrpwv(z)
#println(tfr)
