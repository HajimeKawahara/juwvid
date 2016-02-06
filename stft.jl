module stft

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
    Lh=round(Int64,(hrow-1)/2) ##??

    tfr=zeros(Complex64,N,N) # plane by default
    for icol=1:N
        ti=t[icol]
        # tau=collect(-minimum([round(N/2)-1,Lh,ti-1]):minimum([round(N/2)-1,Lh,xrow-ti]))
        # indices= rem(N+tau,N)+1; 
        # tfr(indices,icol)=x(ti+tau,1).*conj(h(Lh+1+tau));

        taumin=-minimum([floor(N/2)-1,Lh,ti-1])
        taumax=minimum([floor(N/2)-1,Lh,xrow-ti])
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
