module jnufft
#change this path in ccall for libnufft.so

function call_ionufft1d2(xj,fk,iflag=-1,eps=10.0^-24)
#requested frequency xj
#values in the time sequence fk

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

    #(nj,ms,mx,iflag,eps,xj,fk,cj)
    product = ccall((:__ionufft_MOD_ionufft1d2, "/Users/kawahara/juwvid/libnufft.so"),Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Complex128},Ptr{Complex128}),&nj,&ms,&mx,&iflag,&eps,xj0,fk0,cj)
#    return cj[1:ms]
    return unshift!(cj[1:mx-1], cj[mx])
end

function call_ionufft1d3(sk,xj,cj,iflag=-1,eps=10.0^-24)
#requested frequency: sk
#input nonuniform time sequence (time): xj
#values in the time sequence (value): cj

    nk=Int32(length(sk))
    nj=Int32(length(xj))    
    iflag=Int32(iflag)
    eps=Float64(eps)

    fk = zeros(Complex128,nk) 
    sk0 = zeros(Float64,nk)     
    for k = 1:nk
#        sk0[k] = 2*pi*sk[k]/nk
        sk0[k] = sk[k]
    end

    xj0 = zeros(Float64,nj)     
    for j = 1:nj
        xj0[j] = 2*pi*(xj[j]-nj/2-1)/nj
    end

    cj0 = zeros(Complex128,nj) 
    for j = 1:round(Int,nj/2)
        cj0[j] = cj[round(Int,j+nj/2)]
    end
    for j = round(Int,nj/2+1):nj          
        cj0[j] = cj[round(Int,j-nj/2)]
    end

    #ionufft1d3(nj,nk,iflag,eps,xj,sk,cj,fk)
    product = ccall((:__ionufft_MOD_ionufft1d3, "/Users/kawahara/juwvid/libnufft.so"),Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Complex128},Ptr{Complex128}),&nj,&nk,&iflag,&eps,xj0,sk0,cj0,fk)
    
    return unshift!(fk[1:nk-1], fk[nk])

end

function test_1d2(ms=8,nj=16)
##TEST for 1d2
##value for time sequence (discrete position)
println("########### JNUFFT 1D2 ###########")
fk = collect(1:ms)
println("value for time sequence=",fk[1],"-",fk[end])
##arbitrary frequency
xj = collect(1:nj)*real(ms/nj)     
println("requested frequency=",xj[1],"-",xj[end])
cj=jnufft.call_ionufft1d2(xj,fk,-1,10.0^-32)
println("derived coefficient=",cj)
end

function test_1d3(nj=8,nk=16)
#test_1d2()
println("########### JNUFFT 1D3 ###########")
xj= collect(1:nj)
println("time grids=",xj)
cj = collect(1:nj)
println("value for time sequence=",cj[1],"-",cj[end])
##arbitrary frequency
sk = collect(1:nk)*real(nj/nk)     
println("requested frequency=",sk[1],"-",sk[end])
fk=jnufft.call_ionufft1d3(sk,xj,cj,-1,10.0^-24)
println("derived coefficient=",fk)
end


function test_juliafft(nj)
cj = collect(1:nj)
println("########### JULIA FFT ###########")
println("value for time sequence=",cj[1],"-",cj[end])
##arbitrary frequency
fk=fft(cj)
println("derived coefficient=",fk)
end


end
jnufft.test_juliafft(8)
jnufft.test_1d2(8,10)
jnufft.test_1d3(8,10)