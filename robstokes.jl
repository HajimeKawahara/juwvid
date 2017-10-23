module robstokes
import jnufft

function robstpwv(x,y,t=NaN,f=NaN,itc=NaN,silent=0,nwindow=4,Nz=1)
    #method = median : robust Wigner distribution
    xrow = size(x)[1] 
    if isnan.(t)[1] t=collect(1:xrow) end
    N=xrow
    if isnan.(itc)[1]  
        Nt=N
        itc=collect(1:Nt)
    else
        Nt=length(itc)
    end
    hlength=floor(N/nwindow)
    hlength=hlength+1-rem.(hlength,2)        
    Lh=round.(Int,(hlength-1)/2)

    tfrI=zeros(Complex64,N,Nt)
    tfrQ=zeros(Complex64,N,Nt)
    tfrU=zeros(Complex64,N,Nt)
    tfrV=zeros(Complex64,N,Nt)

    for icol=1:Nt
        #ti=t[icol]
        ti=t[itc[icol]]
        taumax=minimum([ti-1,xrow-ti,round.(N/2)-1,Lh])
        tau=round.(Int64,-taumax:taumax); 
        indices=round.(Int64,rem.(N+tau,N)+1)
        tfrI[indices,icol] = (x[ti+tau].*conj(x[ti-tau])+y[ti+tau].*conj(y[ti-tau]))
        tfrQ[indices,icol] = (x[ti+tau].*conj(x[ti-tau])-y[ti+tau].*conj(y[ti-tau]))
        tfrU[indices,icol] = (x[ti+tau].*conj(y[ti-tau])+y[ti+tau].*conj(x[ti-tau]))
        tfrV[indices,icol] = -(y[ti+tau].*conj(x[ti-tau])-x[ti+tau].*conj(y[ti-tau]))*im

        tau=round.(N/2);
        if ti<=xrow-tau && ti>=tau+1 && tau<=Lh
            tfrI[tau+1,icol] = 0.5*(x[ti+tau]*conj(x[ti-tau]) + x[ti-tau]*conj(x[ti+tau])) + 0.5*(y[ti+tau]*conj(y[ti-tau]) + y[ti-tau]*conj(y[ti+tau]))
            tfrQ[tau+1,icol] = 0.5*(x[ti+tau]*conj(x[ti-tau]) + x[ti-tau]*conj(x[ti+tau])) - 0.5*(y[ti+tau]*conj(y[ti-tau]) + y[ti-tau]*conj(y[ti+tau]))
            tfrU[tau+1,icol] = 0.5*(x[ti+tau]*conj(y[ti-tau]) + x[ti-tau]*conj(y[ti+tau])) + 0.5*(y[ti+tau]*conj(x[ti-tau]) + y[ti-tau]*conj(x[ti+tau]))
            tfrV[tau+1,icol] = -(0.5*(y[ti+tau]*conj(x[ti-tau]) + y[ti-tau]*conj(x[ti+tau])) - 0.5*(x[ti+tau]*conj(y[ti-tau]) + x[ti-tau]*conj(y[ti+tau])))*im
        end
    end

    if silent==0 println("Robust Stokes distribution (very slow)") end
    if isnan.(f)[1] 
        f=collect(linspace(1,Nt,Nt))
    end
    Nf=size(f)[1]
    m=collect(1:Nt)
    tfrnewI=zeros(Complex64,Nf,Nt) 
    tfrnewQ=zeros(Complex64,Nf,Nt) 
    tfrnewU=zeros(Complex64,Nf,Nt) 
    tfrnewV=zeros(Complex64,Nf,Nt) 

    for i=1:Nt
        if rem.(i,256)==1 println(i,"/",Nt) end
        for j=1:Nf
            arr=real(tfrI[:,i].*exp.(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            tfrnewI[j,i]=median(arr[arr.!=0.0])
            
            arr=real(tfrQ[:,i].*exp.(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            tfrnewQ[j,i]=median(arr[arr.!=0.0])

            arr=real(tfrU[:,i].*exp.(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            tfrnewU[j,i]=median(arr[arr.!=0.0])

            arr=real(tfrV[:,i].*exp.(-2.0*pi*im*(m[:]-1)*(f[j]-1)/N))
            tfrnewV[j,i]=median(arr[arr.!=0.0])

        end            
    end
    return tfrnewI, tfrnewQ, tfrnewU, tfrnewV 
end

end

#import DSP
#y=linspace(0.0,16.0,16)
#z=DSP.Util.hilbert(y)
#tfr=cohenclass.tfrwv(z)
#nsample=16
#tfr=cohenclass.tfrwv(z,NaN,NaN,[1,2],NaN,0)
##tfr=cohenclass.tfrpwv(z)
#println(tfr)

