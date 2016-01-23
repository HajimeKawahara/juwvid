module pwvaw #pseudo Wigner distribution with the adaptive algorithm
import jnufft

function awpwv(x,varrat,y=NaN,t=NaN,f=NaN,N=NaN,kappa=5.0,delkappa=0.97,silent=0,method="mean",use_nufft=true)
    # adaptive algorithm of pseudo Wigner-Ville distribution 
    # Reference: Time Frequency Signal Analysis with Applications, 2013, Stanlovi, Dakovic, and Thayaparan p420
    # varrat[1:length(x)]: variance ratio of the noise to the signal (sigma_e/A)^2
    # If f is not NaN, use NUFFT by default

    wtype="rect"
#    wtype="Hamming"

    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end
    if isnan(N) N=xrow end
    if isnan(y)[1] 
        if silent ==0 println("Single pseudo Wigner Ville with the Adaptive Algorithm") end
        y=x
        sw=0
    else
        if silent ==0 println("Cross pseudo Wigner Ville with the Adaptive Algorithm") end
        sw=1
    end

    #prepare tfrnew
    if isnan(f)[1]
        Nf=N
    else
        Nf=size(f)[1]    
    end
    tfrnew=zeros(Complex64,Nf,N) 

    #set window array
    M=N/2
    harray=[M]
    while M>4
        M=M/2
        push!(harray,M)       
    end
    harray=harray[end:-1:1]
    windows=[]

    #
    for icol=1:N
        crit=true
        prevcrit=false
        ih=0

        while crit && ih < length(harray) 
            ih=ih+1
            tfrvec=zeros(Complex64,N) #
            hlength=floor(harray[ih])
            hlength=hlength+1-rem(hlength,2)        
            
            if ismatch(r"Hamming",wtype)         
                h=0.54 - 0.46*cos(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
            else
                h=ones(round(Int,hlength)) #Rectangular
            end    

            hrow = size(h)[1] 
            Lh=round(Int,(hrow-1)/2)
            h=h/h[Lh+1] # normalization so that h(0) = 1.0
            ti=t[icol]
            taumax=minimum([ti-1,xrow-ti,round(N/2)-1,Lh])
            tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
            tfrvec[indices] = h[Lh+1+tau].*x[ti+tau].*conj(y[ti-tau])
            tau=round(N/2); 
            if ti<=xrow-tau && ti>=tau+1 && tau<=Lh
                tfrvec[tau+1] = 0.5*(h[Lh+1+tau]*x[ti+tau]*conj(y[ti-tau]) + h[Lh+1-tau]*x[ti-tau]*conj(y[ti+tau]))
            end

            if isnan(f)[1] 
                tfrvecp=fft(tfrvec) #fft
            elseif use_nufft
                tfrvecp=jnufft.call_ionufft1d2(f,tfrvec,-1,10.0^-28)[1:Nf] #nufft
            else
                m=collect(1:N)
                for j=1:Nf                
                    tfrvecp[j]=sum(tfrvec.*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N)) #dft
                end            
            end
            fhs=indmax(abs(tfrvecp))            
            sigmahs=sqrt((kappa+delkappa)*varrat[icol]/(hlength^3))

            if ih>1 && abs(fhsprev - fhs)/(2.0*N) > (kappa+delkappa)*(sigmahsprev+sigmahs)      
                tfrnew[:,icol]=tfrvecprev
                push!(windows,harray[ih-1])
                crit=false
            elseif ih==length(harray)
                tfrnew[:,icol]=tfrvecprev
                push!(windows,harray[ih-1])
                crit=false
            end            

            fhsprev=fhs
            sigmahsprev=sigmahs
            tfrvecprev=copy(tfrvecp)
        end        
    end
    return tfrnew, windows
end

end

#import DSP
#import cohenclass
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=cohenclass.tfrpwv(z)
#tfr,windows=pwvaw.awpwv(z)
