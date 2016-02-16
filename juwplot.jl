module juwplot
import PyPlot

function tfrshow(tfrs,dx,x1,xend,fin1=NaN,finend=NaN,asp=0.7,cmap="CMRmap",vmin=NaN,vmax=NaN)
    nsample=size(tfrs,2)
    freqfac=1/nsample/dx/2
    if isnan(fin1) || isnan(finend)
        f1=freqfac
        fe=freqfac*size(tfrs,1)
    else
        offset=(finend-fin1)/nsample
        f1=(fin1-offset)*freqfac
        fe=(finend-offset)*freqfac            
    end

    if isnan(vmin) || isnan(vmax)
        a=PyPlot.imshow(real(tfrs[end:-1:2,:]),cmap=cmap,extent=[x1,xend,f1,fe],aspect=asp*(xend-x1)/(fe-f1),interpolation="nearest")
    else
        a=PyPlot.imshow(real(tfrs[end:-1:2,:]),cmap=cmap,extent=[x1,xend,f1,fe],aspect=asp*(xend-x1)/(fe-f1),interpolation="nearest",vmin=vmin,vmax=vmax)
    end

    return a
end

function wtfrshow(tfrs,dx,x1,xend,fin1=NaN,finend=NaN,asp=0.7,cmap="CMRmap",vmin=NaN,vmax=NaN)
    nsample=size(tfrs,2)
    freqfac=1/nsample/dx/2
    if isnan(fin1) || isnan(finend)
        f1=freqfac
        fe=freqfac*size(tfrs,1)
    else
        offset=(finend-fin1)/nsample
        f1=(fin1-offset)*freqfac
        fe=(finend-offset)*freqfac            
    end
    endi=round(Int,size(tfrs,1)/2.0)

    if isnan(vmin) || isnan(vmax)
        a=PyPlot.imshow(real(tfrs[endi:-1:2,:]),cmap=cmap,extent=[x1,xend,f1,fe],aspect=asp*(xend-x1)/(fe-f1),interpolation="nearest")
    else
        a=PyPlot.imshow(real(tfrs[endi:-1:2,:]),cmap=cmap,extent=[x1,xend,f1,fe],aspect=asp*(xend-x1)/(fe-f1),interpolation="nearest",vmin=vmin,vmax=vmax)
    end

    return a
end

end