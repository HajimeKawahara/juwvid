module stokes
using stft
using cohenclass
using smethod
import DSP

function get_StokesTFR(W11,W22,W12,W21)
    I=W11+W22
    Q=W11-W22
    U=W12+W21
    V=(W21-W12)*(im)
    return I,Q,U,V
end

function get_StokesTFR_LR(W11,W22,W12,W21)
    I=W11+W22
    V=W11-W22
    Q=W12+W21
    U=(W21-W12)*(im)
    return I,Q,U,V
end



function wvstokes(y1,y2,t=NaN,f=NaN,itc=NaN,silent=0,method="mean",def="wolf")
    #Stokes Distribution for Wigner Ville Distribution
    z1=DSP.Util.hilbert(y1)
    z2=DSP.Util.hilbert(y2)

    W11=cohenclass.tfrwv(z1,NaN,t,f,itc,silent,method)
    W22=cohenclass.tfrwv(z2,NaN,t,f,itc,silent,method)
    W12=cohenclass.tfrwv(z2,z1,t,f,itc,silent,method)
    W21=cohenclass.tfrwv(z1,z2,t,f,itc,silent,method)

    if def == "wolf"
        I,Q,U,V=get_StokesTFR(W11,W22,W12,W21)
    elseif def == "LR"
        I,Q,U,V=get_StokesTFR_LR(W11,W22,W12,W21)
    end

    return I,Q,U,V
end

function pwvstokes(y1,y2,t=NaN,f=NaN,itc=NaN,h=NaN,silent=0,method="fft",nwindow=4,Nz=1,def="wolf")
    #Stokes Distribution for pseudo Wigner Ville Distribution
    if def=="wolf"
        z1=DSP.Util.hilbert(y1)
        z2=DSP.Util.hilbert(y2)
    elseif def == "LR" 
        z1=(y1)
        z2=(y2)
    end

    W11=cohenclass.tfrpwv(z1,NaN,t,f,itc,h,silent,method,nwindow,Nz)
    W22=cohenclass.tfrpwv(z2,NaN,t,f,itc,h,silent,method,nwindow,Nz)
    W12=cohenclass.tfrpwv(z2,z1,t,f,itc,h,silent,method,nwindow,Nz)
    W21=cohenclass.tfrpwv(z1,z2,t,f,itc,h,silent,method,nwindow,Nz)

    if def == "wolf"
        I,Q,U,V=get_StokesTFR(W11,W22,W12,W21)
    elseif def == "LR" 
        I,Q,U,V=get_StokesTFR_LR(W11,W22,W12,W21)
    end

    return I,Q,U,V
end

function smstokes(y1,y2,Lp=6,f=NaN,silent=0,nwindow=4,fps=NaN,fpe=NaN,itc=NaN,def="wolf")
    #Stokes Distribution for S-method
    sm11=smethod.tfrsm(y1,NaN,Lp,f,nwindow,silent,fps,fpe,itc)
    sm22=smethod.tfrsm(y2,NaN,Lp,f,nwindow,silent,fps,fpe,itc)
    sm12=smethod.tfrsm(y2,y1,Lp,f,nwindow,silent,fps,fpe,itc)
    sm21=smethod.tfrsm(y1,y2,Lp,f,nwindow,silent,fps,fpe,itc)
    if def == "wolf"
        I,Q,U,V=get_StokesTFR(sm11,sm22,sm12,sm21)
    elseif def == "LR"
        I,Q,U,V=get_StokesTFR_LR(sm11,sm22,sm12,sm21)
    end

    return I,Q,U,V
end

end
