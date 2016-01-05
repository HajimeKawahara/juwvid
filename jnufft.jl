module jnufft
#change the path in ccall for libnufft.so

function call_ionufft1d2(xj,fk,iflag=-1,eps=10.0^-24)

    nj=Int32(length(xj))    
#    println("nj=",nj)
    ms=Int32(length(fk))
#    println("ms=",ms)
    iflag=Int32(iflag)
    eps=Float64(eps)

    mx=maximum([nj,ms])
    cj = zeros(Complex128,mx) 

    xj0 = zeros(Float64,nj)     
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
    product = ccall((:__ionufft_MOD_ionufft1d2, "/Users/kawahara/juwvid/libnufft.so"),Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Complex128},Ptr{Complex128}),&nj,&ms,&mx,&iflag,&eps,xj0,fk0,cj)
#    return cj[1:ms]
    return unshift!(cj[1:mx-1], cj[mx])
end

end

##value for discrete position
#ms=8
#fk = collect(1:ms)
##arbitrary frequency
#nj=120
#xj = collect(1:nj)*real(ms/nj)     
#cj=jnufft.call_ionufft1d2(xj,fk,-1,10.0^-32)
#println(cj)

