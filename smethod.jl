module smethod
using stft 


function tfrsm(x,y=NaN,Lp=6,f=NaN,nwindow=4,silent=0,fps=NaN,fpe=NaN,itc=NaN)
    #
    #normal/cross S-method
    #fps: pick-up frequency (start)
    #fpe: pick-up frequency (end)
    if isnan.(y)[1]
        if silent ==0  println("Single S-method") end
        sw=0
    else 
        if silent==0 println("Cross S-method") end
        sw=1
    end

    if isnan.(itc)[1]  
        itc=collect(1:length(x))
    end

    if isnan.(f)[1]  
        tfrstftx=stft.tfrstft(x,NaN,NaN,NaN,itc,NaN,nwindow)    
        if sw==0
            tfrstfty=tfrstftx
        else
            tfrstfty=stft.tfrstft(y,NaN,NaN,NaN,itc,NaN,nwindow)    
        end
    else
        tfrstftx=stft.tfrstft(x,NaN,NaN,f,itc,NaN,nwindow)            
        if sw==0
            tfrstfty=tfrstftx
        else
            tfrstfty=stft.tfrstft(y,NaN,NaN,f,itc,NaN,nwindow)            
        end
    end

    nsamplef=size(tfrstftx)[1]
    nsamplet=size(tfrstftx)[2]

    if isnan.(fps) || isnan.(fpe) 
        ks=1
        ke=nsamplef 
        sm=zeros(Complex64,nsamplef,nsamplet)
    else
        ks=maximum([1,2*fps-Lp])
        ke=minimum([nsamplef,2*fpe+Lp])
        sm=zeros(Complex64,ke-ks+1,nsamplet)
    end

#    for k=1:nsamplef
    for k=ks:ke
        kq=k-ks+1
        sm[kq,:]=(tfrstftx[k,:].*conj(tfrstfty[k,:]))
        for i=1:Lp
            if k+i<=nsamplef && k-i>0
                sm[kq,:]=sm[kq,:]+(tfrstftx[k+i,:].*conj(tfrstfty[k-i,:])+tfrstftx[k-i,:].*conj(tfrstfty[k+i,:]))
            end
        end
    end

    return sm
end


function tfrsm_single(x,Lp,f=NaN,nwindow=4,fps=NaN,fpe=NaN,itc=NaN)
    #only for single S-method (but, a bit faster)
    #fps: pick-up frequency (start)
    #fpe: pick-up frequency (end)
    if isnan.(itc)[1]  
        itc=collect(1:length(x))
    end

    if isnan.(f)[1]  
        tfrstft=stft.tfrstft(x,NaN,NaN,NaN,itc,NaN,nwindow)    
    else
        tfrstft=stft.tfrstft(x,NaN,NaN,f,itc,NaN,nwindow)            
    end

    nsamplef=size(tfrstft)[1]
    nsamplet=size(tfrstft)[2]

    if isnan.(fps) || isnan.(fpe) 
        ks=1
        ke=nsamplef 
        sm=zeros(nsamplef,nsamplet)
    else
        ks=maximum([1,2*fps-Lp])
        ke=minimum([nsamplef,2*fpe+Lp])
        sm=zeros(ke-ks+1,nsamplet)
    end

#    for k=1:nsamplef
    for k=ks:ke
        kq=k-ks+1
        sm[kq,:]=(tfrstft[k,:].*conj(tfrstft[k,:]))
        for i=1:Lp
            if k+i<=nsamplef && k-i>0
                sm[kq,:]=sm[kq,:]+2*real(tfrstft[k+i,:].*conj(tfrstft[k-i,:]))
            end
        end
    end

    return sm
end


end
