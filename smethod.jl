module smethod
using stft 

function tfrsm(x,Lp,f=NaN,nwindow=4)

    if isnan(f)[1] 
        tfrstft=stft.tfrstft(x,NaN,NaN,NaN,NaN,nwindow)    
    else
        tfrstft=stft.tfrstft(x,NaN,NaN,f,NaN,nwindow)            
    end

    nsamplef=size(tfrstft)[1]
    nsamplet=size(tfrstft)[2]
    sm=zeros(nsamplef,nsamplet)

    for k=1:nsamplef
        sm[k,:]=(tfrstft[k,:].*conj(tfrstft[k,:]))
        for i=1:Lp
            if k+i<=nsamplef && k-i>0
                sm[k,:]=sm[k,:]+2*real(tfrstft[k+i,:].*conj(tfrstft[k-i,:]))
            end
        end
    end

    return sm
end

end