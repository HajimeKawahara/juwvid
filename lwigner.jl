module lwigner
using smethod

function trflw2L(tfr,Lp)

    nsamplef=size(tfr)[1]
    nsamplet=size(tfr)[2]
    trflw=zeros(nsamplef,nsamplet)

    for k=1:nsamplef
        trflw[k,:]=(tfr[k,:].*tfr[k,:])
        for i=1:Lp
            if k+i<=nsamplef && k-i>0
                trflw[k,:]=trflw[k,:]+2*(tfr[k+i,:].*(tfr[k-i,:]))
            end
        end
    end

    return trflw

end

end