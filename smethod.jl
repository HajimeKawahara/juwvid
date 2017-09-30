module smethod
using stft 

function tfrsm(x,Lp,f=NaN,nwindow=4,fps=NaN,fpe=NaN,itc=NaN)

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
