module pwvaw #pseudo Wigner distribution with the adaptive 

function awpwv(x,y=NaN,t=NaN,f=NaN,N=NaN,h=NaN,silent=0,method="mean",nwindow=4,use_nufft=true)
    #method = median : robust Wigner distribution

    wtype="rect"

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

    #prepare tfrnew
    if isnan(f)
        Nf=N
    else
        Nf=size(f)[1]    
    end
    tfrnew=zeros(Complex64,Nf,N) 

    #set window array
    M=N/4
    harray=[M]
#    println(harray)
    while M>4
        M=M/2
        push!(harray,M)       
    end
    harray=harray[end:-1:1]
    windows=[]
#    println("Adaptive window from =",harray)

    for icol=1:N
        ih=0
        crit=true
        prevcrit=false
        while crit && ih < length(harray) 
            ih=ih+1
#            println("ih=",ih)
            tfrvec=zeros(Complex64,N) #
            hlength=floar(harray[ih])
            hlength=hlength+1-rem(hlength,2)        
            
#            if ismatch(r"Hamming",wtype)         
                h=0.54 - 0.46*cos(2.0*pi*(1:hlength)/(hlength+1)) #Hamming
#            else
#                #rectangular
#                h=ones(round(Int,hlength))
#            end    

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
#                if silent==0 println("Use fft.") end
                tfrnew[:,icol]=fft(tfrvec)
            elseif use_nufft
#                if silent==0 println("Use nufft.") end
                tfrnew[:,icol]=jnufft.call_ionufft1d2(f,tfrvec,-1,10.0^-28)[1:Nf]
            else
#                if silent==0 println("Use dfft.") end 
                m=collect(1:N)
                for j=1:Nf                
                    tfrnew[j,icol]=sum(tfrvec.*exp(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
                end            
            end

            fhs=indmax(abs(tfrvec))            
            ##### CRITERION #####
            A=1.0
            sigmae=0.00
            kappa=5
            delkappa=0.97
            p=0.72
            p1=0.97
            #####################
            sigmahs=sqrt(6.0*sigmae^2/(A^2*hlength^3))
            if ih==5                
#            if ih>1 && abs(fhsprev - fhs) <= (kappa+delkappa)*(sigmahsprev+sigmahs)                
                #println("change ","ih=",ih,"left:",abs(fhsprev - fhs),"right:",(kappa+delkappa)*(sigmahsprev+sigmahs))
                prevcrit=true
                if ih == length(harray)
                    crit=false
                    push!(windows,harray[ih])
                end

            end
 
            if prevcrit
                crit=false
                push!(windows,harray[ih])
            end
            
            fhsprev=fhs
            sigmahsprev=sigmahs

        end
        
    end
    println(windows)
    return tfrnew, windows

end

end

#import DSP
#nsample=1024
#x=linspace(0.0,nsample,nsample)
#y=x
#z=DSP.Util.hilbert(y)
#tfr,windows=pwvaw.awpwv(z)
#println(windows)

