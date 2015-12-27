module polywv

function tfr4powv(x,t=NaN,N=NaN)
    #fourth order polinomial
    xcol,xrow = size(x) #xcol for cross-rwv
    if isnan(t) t=(1:xrow)' end
    if isnan(N) N=xrow end
    #information
    println("Sizes of x and t")
    println(size(x))
    println(size(t))

    #error handling
    if N<0; println("N must be greater than zero"); exit(); end
    if xcol==0 || xcol>2; println("X must have one or two columns"); exit() end
    if xcol==1; println("Single Wigner Ville"); end
    if xcol==2; println("Cross Wigner Ville"); end
    if nextpow2(N)!=N; println("For a faster computation, N should be a power of two\n"); end

    tfr=zeros(Complex64,N,N) # plane by default

    for icol=1:N
        ti=t[icol]
        taumax=minimum([ti-1,xrow-ti,round(N/4)-1])
        tau=round(Int64,-taumax:taumax); indices=round(Int64,rem(N+tau,N)+1)
        tfr[indices,icol] = x[1,ti+tau].*conj(x[xcol,ti-tau]).*x[1,ti+tau].*conj(x[xcol,ti-tau]) #                

        tau=round(N/4); 
        if ti<=xrow-tau && ti>=tau+1
            tfr[tau+1,icol] = 0.5*(x[1,ti+tau]^2*conj(x[xcol,ti-tau])^2 + x[1,ti-tau]^2*conj(x[xcol,ti+tau]^2))
        end
    end
    for i=1:N
        tfr[:,i]=fft(tfr[:,i])
    end
    if xcol==1
        tfr=real(tfr)
    end

    return tfr
end

end

import DSP
y=[1.,2.,3.,5.,1.,2.,3.,5.,1.,2.,3.,5.,1.,2.,3.,5.]
ya=DSP.Util.hilbert(y) # transpose is necessary 
y=conj(ya')
tfr=tfr4powv(y)
println(real(tfr))
