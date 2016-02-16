module juwplot
import PyPlot

function tfrshow(tfrs,dx,x1,xend,fin1=NaN,finend=NaN,asp=0.7,cmap="CMRmap")
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
    a=PyPlot.imshow(real(tfrs[end:-1:2,:]),cmap=cmap,extent=[x1,xend,f1,fe],aspect=asp*(xend-x1)/(fe-f1),interpolation="nearest")
    return a
end

end
