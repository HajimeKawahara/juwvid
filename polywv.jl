module polywv
import jnufft
import stft
import smethod
import lwigner

function tfrpowv(x,y=NaN,t=NaN,f=NaN,nwindow=2,Lp=NaN,Lpsm=4,Lpl2=4,Lppwv=4,silent=0)
    # Boashash 15 eq(6.2.19) p349
    xrow = size(x)[1] 
    if isnan(t)[1] t=collect(1:xrow) end

    if isnan(Lp) 
        if isnan(Lppwv) || isnan(Lpsm) || isnan(Lpl2)
            println("Missing Lp.")
            return NaN
        end
        Lp=Lppwv
    else
        Lppwv=Lp
        Lpsm=Lp
        Lpl2=Lp
    end

    nsamplet=length(x)
    #alias free sm
    if isnan(f)[1]
        afwv=smethod.tfrsm(x,Lpsm,NaN,nwindow)
        nsamplef=nsamplet
    else
        afwv=smethod.tfrsm(x,Lpsm,f,nwindow)
        nsamplef=length(f)
    end

    #L2 wigner
    tfrlw2=lwigner.tfrlw2L(afwv,Lpl2)       

    A=0.85/1.35
    tfrpwv=zeros(nsamplef,nsamplet)

    for k=1:nsamplef
        for i=-Lp:Lp
            j=round(Int,float(i)/A)
            if k+j<=nsamplef && k+j >0                
                tfrpwv[k,:]=tfrpwv[k,:]+tfrlw2[k+i,:].*afwv[k+j,:]
            end
        end
    end
    
    return tfrpwv
    
end
end

#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y) # transpose is necessary 
#tfr=polywv.tfrpowv(z)
