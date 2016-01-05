module jnufft

function call_ionufft1d2(xj,fk,iflag=-1,eps=10.0^-24)

    nj=Int32(length(fk))
    ms=Int32(length(xj))    
    iflag=Int32(iflag)
    eps=Float64(eps)

    cj = zeros(Complex128,nj) 

    xj0 = zeros(Float64,ms)     
    for j = 1:nj
        xj0[j] = 2*pi*xj[j]/ms
    end

    fk0 = zeros(Complex128,ms) 
    for j = 1:round(Int,ms/2)
        fk0[j] = fk[round(Int,j+ms/2)]
    end
    for j = round(Int,ms/2+1):ms          
        fk0[j] = fk[round(Int,j-ms/2)]
    end
    product = ccall((:__ionufft_MOD_ionufft1d2, "libnufft.so"),Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Complex128},Ptr{Complex128}),&nj,&ms,&iflag,&eps,xj0,fk0,cj)
    #    return cj[1:ms]
    return unshift!(cj[1:ms-1], cj[ms])
end

end

#ms=8
#nj=8 
#xj = collect(1:ms)     
#fk = collect(1:ms)
#cj=jnufft.call_ionufft1d2(xj,fk,-1,10.0^-32)
#println(cj)

