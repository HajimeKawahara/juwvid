module smethod
using stft 

function tfrsm(x,Lp,nwindow=4)
    nsample=length(x)
    tfrstft=stft.tfrstft(x,NaN,NaN,NaN,nwindow)    
    sm=zeros(nsample,nsample);
    for k=1:nsample
        sm[k,:]=(tfrstft[k,:].*conj(tfrstft[k,:]))
        for i=1:Lp
            if k+i<=nsample && k-i>0
                sm[k,:]=sm[k,:]+2*real(tfrstft[k+i,:].*conj(tfrstft[k-i,:]))
            end
        end
    end

    return sm
end

end