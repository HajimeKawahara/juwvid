module estif
#import DSP
import cohenclass
import extif
import polywv

function ifestxvwd(z,ynorm,dx,niter)
    tfrn=cohenclass.tfrwv(z,NaN,NaN,NaN,1)
    indfn=extif.maxif(real(tfrn))
    for it=1:niter        
        zhat=exp(im*cumsum(indfn)*ynorm*dx)
        tfrn=cohenclass.tfrwv(z,zhat,NaN,NaN,1)
        indfn=extif.maxif(abs(tfrn))
    end
    return indfn
end

function ifestxpvwd(z,ynorm,dx,niter)
    tfrn=cohenclass.tfrpwv(z,NaN,NaN,NaN,NaN,1)
    indfn=extif.maxif(real(tfrn))
    for it=1:niter
        zhat=exp(im*cumsum(indfn*ynorm*dx))
        tfrn=cohenclass.tfrpwv(z,zhat,NaN,NaN,NaN,1)
        indfn=extif.maxif(abs(tfrn))
    end
    return indfn
end

function ifestxpowv(z,ynorm,dx,niter, tfrsw=false)
    tfrn=polywv.tfrpowv(z,NaN,NaN,NaN,1)
    indfn=extif.maxif(real(tfrn))
    for it=1:niter        
        zhat=exp(im*cumsum(indfn)*ynorm*dx)
        tfrn=polywv.tfrpowv(z,zhat,NaN,NaN,1)
        indfn=extif.maxif(abs(tfrn))
    end

    if tfrsw
        return indfn, tfrn
    else
        return indfn
    end
end

end