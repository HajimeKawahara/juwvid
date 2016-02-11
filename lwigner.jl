module lwigner
using smethod

function tfrlw2L(tfr,Lplw)

    nsamplef=size(tfr)[1]
    nsamplet=size(tfr)[2]
    tfrlw=zeros(nsamplef,nsamplet)

    for k=1:nsamplef
        tfrlw[k,:]=(tfr[k,:].*tfr[k,:])
        for i=1:Lplw
            if k+i<=nsamplef && k-i>0
                tfrlw[k,:]=tfrlw[k,:]+2*(tfr[k+i,:].*(tfr[k-i,:]))
            end
        end
    end

    return tfrlw

end

end