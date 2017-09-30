module cohenclass
import jnufft

function tfrwv(x,y=NaN,t=NaN,f=NaN,itc=NaN,silent=0,method="mean",use_nufft=true)
    xrow = size(x)[1] 

    if isnan.(t)[1] t=collect(1:xrow) end
    N=xrow
    if isnan.(itc)[1]  
        Nt=N
        itc=collect(1:Nt)
    else
        Nt=length(itc)
    end

    if isnan.(y)[1]
        if silent ==0  println("Single Wigner Ville") end
        y=x
        sw=0
    else 
        if silent==0 println("Cross Wigner Ville") end
        sw=1
    end
    tfr=zeros(Complex64,N,Nt) # plane by default
    for icol=1:Nt       
#        ti=t[icol]
        ti=round.(Int64,t[itc[icol]])
        taumax=minimum([ti-1,xrow-ti,round.(N/2)-1])
        tau=round.(Int64,-taumax:taumax); indices=round.(Int64,rem.(N+tau,N)+1)
        tfr[indices,icol] = x[ti+tau].*conj(y[ti-tau]) #                
        tau=round.(Int64,N/2); 
        if ti<=xrow-tau && ti>=tau+1 
            tfr[tau+1,icol] = 0.5*(x[ti+tau]*conj(y[ti-tau]) + x[ti-tau]*conj(y[ti+tau]))
        end
    end


    ###Choose FFT or DFT
    if isnan.(f)[1] && ismatch(r"mean",method)
        if silent==0 println("Use fft.") end
        for i=1:Nt
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    elseif ismatch(r"mean",method)
        if silent==0 println("Use nufft.") end
        Nf=size(f)[1]
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
            tfrnew[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-28)[1:Nf]
        end        
        return tfrnew
	
    elseif ismatch(r"nufft",method)
        if silent==0 println("Use nufft.") end
        Nf=size(f)[1]
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
            tfrnew[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-28)[1:Nf]
        end        
        return tfrnew

    elseif isnan.(f)[1] && ismatch(r"fft",method)
        if silent==0 println("Use fft.") end
        for i=1:Nt
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr

    elseif ismatch(r"dft",method)
        if silent==0 println("Use Direct DFT.") end
        Nf=size(f)[1]
        m=collect(1:Nt)
        tfrnew=zeros(Complex64,Nf,Nt)        
        for i=1:Nt
            for j=1:Nf                
                tfrnew[j,i]=sum(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            end            
        end
        return tfrnew
    elseif ismatch(r"robust",method)
        if silent==0 println("Robust Wigner distribution") end
        Nf=size(f)[1]
        m=collect(1:Nt)
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
            if rem.(i,256)==1 println(i,"/",Nt) end
            for j=1:Nf             
                arr=real(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
                arr=arr[arr.>0.0]
                tfrnew[j,i]=median(arr)
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

function tfrpwv(x,y=NaN,t=NaN,f=NaN,itc=NaN,h=NaN,silent=0,method="fft",nwindow=4,Nz=1)
    #method = median : robust Wigner distribution
    xrow = size(x)[1] 
    if isnan.(t)[1] t=collect(1:xrow) end
    N=xrow
    if isnan.(itc)[1]  
        Nt=N
        itc=collect(1:Nt)
    else
        Nt=length(itc)
    end

    if isnan.(y)[1] 
        if silent ==0 println("Single pseudo Wigner Ville") end
        y=x
        sw=0
    else
        if silent ==0 println("Cross pseudo Wigner Ville") end
        sw=1
    end
    if isnan.(h)[1] 
        hlength=floor(N/nwindow)
        hlength=hlength+1-rem.(hlength,2)        
        h=0.54 - 0.46*cos.(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
    end    
    hrow = size(h)[1] 
    Lh=round.(Int,(hrow-1)/2)
    h=h/h[Lh+1]

    tfr=zeros(Complex64,N,Nt) # plane by default

    for icol=1:Nt
        #ti=t[icol]
        ti=t[itc[icol]]
        taumax=minimum([ti-1,xrow-ti,round.(N/2)-1,Lh])
        tau=round.(Int64,-taumax:taumax); 
        indices=round.(Int64,rem.(N+tau,N)+1)
        tfr[indices,icol] = h[Lh+1+tau].*x[ti+tau].*conj(y[ti-tau])
        tau=round.(N/2); 
        if ti<=xrow-tau && ti>=tau+1 && tau<=Lh
            tfr[tau+1,icol] = 0.5*(h[Lh+1+tau]*x[ti+tau]*conj(y[ti-tau]) + h[Lh+1-tau]*x[ti-tau]*conj(y[ti+tau]))
        end
    end

    if isnan.(f)[1] && ismatch(r"mean",method)
        #old 
        if silent==0 println("Use fft. (this is the old mode. Use method=fft instead.)") end
        for i=1:Nt
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    elseif ismatch(r"mean",method)
        #old
        if silent==0 println("Use nufft. (this is the old mode. Use method=nufft instead.)") end
        Nf=size(f)[1]
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
                tfrnew[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-28)[1:Nf]
        end
        return tfrnew
    elseif ismatch(r"nufft",method)
        if silent==0 println("Use nufft.") end
        Nf=size(f)[1]
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
                tfrnew[:,i]=jnufft.call_ionufft1d2(f,tfr[:,i],-1,10.0^-28)[1:Nf]
        end
        return tfrnew
    elseif ismatch(r"fft",method) 
        if silent==0 println("Use fft.") end
        for i=1:Nt
            tfr[:,i]=fft(tfr[:,i])
        end
        return tfr
    elseif ismatch(r"zeropad",method) 
        if silent==0 println("Use Zero-padding fft.") end        
        if isnan.(f)[1] 
            Nf=Nt*Nz
            ifs=1
            ife=Nf
        else
            ifs=round.(Int,f[1]*Nz)
            ife=round.(Int,f[end]*Nz)
            Nf=(ife-ifs+1) 
        end
        tfrnew=zeros(Complex64,Nf,Nt) 
        Nh=round.(Int,Nt/2)
        for i=1:Nt
            paddata=vcat(tfr[1:Nh-1,i],zeros(Complex64,(Nz-1)*Nt))
            paddata=vcat(paddata,tfr[Nh:Nt,i])
            tfrnew[:,i]=fft(paddata)[ifs:ife]
        end
        return tfrnew
    elseif ismatch(r"dft",method)
        if silent==0 println("Use dft.") end
        Nf=size(f)[1]
        m=collect(1:Nt)
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
            for j=1:Nf                
                tfrnew[j,i]=sum(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            end            
        end
    elseif ismatch(r"robust",method)
        if silent==0 println("Robust pseudo Wigner distribution") end
        Nf=size(f)[1]
        m=collect(1:Nt)
        tfrnew=zeros(Complex64,Nf,Nt) 
        for i=1:Nt
            if rem.(i,256)==1 println(i,"/",Nt) end
            for j=1:Nf             
                arr=real(tfr[:,i].*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
                arr=arr[arr.>0.0]
                tfrnew[j,i]=median(arr)
            end            
        end
        return tfrnew
    else
        if silent==0 println("No method (fft, zeropad, nufft, dft, robust) exists.") end
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

