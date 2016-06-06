module pm
import juwutils

function compute_all_dphi(x,phase)
    #computing the phase delay matrix
    dx=x[2]-x[1]
    nsample=size(phase)[1]
    offphase=zeros((nsample,nsample))
    freqarr=juwutils.index_to_frequency(collect(1:nsample),NaN,dx,nsample)
    for ifreq=1:nsample
        for it=1:nsample
            offphase[it,ifreq]=mod(2*pi*freqarr[ifreq]*x[it]+pi,2*pi)-pi
        end
    end
    dphiarr=mod(phase[:,:]-offphase[:,:]+pi,2*pi)-pi;

    return dphiarr
end

function compute_dphi(selind,x,phase)
# compute the phase delay
    dx=x[2]-x[1]
    nsample=size(phase)[1]
    selfreq=juwutils.index_to_frequency([selind],NaN,dx,nsample)[1];
    offp=mod(2*pi*selfreq*x[:].+pi,2*pi).-pi
    dphi=mod(transpose(phase[selind,:]).-offp.+pi,2*pi).-pi;
    return dphi
end

end