#simple julia wrapper of fortran 90 subroutine printint in ionufft module 

function call_ionufft1d2(nj,ms,iflag,eps,mx,xj,fk)
#    integer, intent(in) :: ms,nj,iflag
#    integer, intent(in) :: mx !=10000
#    complex*16,intent(in) :: fk(mx)
#    real*8,intent(in) :: xj(mx)
#    real*8,intent(in) :: eps
#    complex*16,intent(out) :: cj(mx)
#    integer :: ier
#    (nj,ms,iflag,eps,mx,xj,fk,cj)

    nj=Int32(nj)
    ms=Int32(ms)
    iflag=Int32(iflag)
    eps=Float64(eps)
    mx=Int32(mx)
    cj = zeros(Complex128,mx) 

    product = ccall((:__ionufft_MOD_ionufft1d2, "nufft_src/libnufft.so"),Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Complex128},Ptr{Complex128}),&nj,&ms,&iflag,&eps,&mx,xj,fk,cj)
#    product = ccall((:__ionufft_MOD_ionufft1d2, "nufft_src/libnufft.so"),Int32,(Ptr{Int32},),nj)

  return cj[1:ms]

end

ms=8
nj=8 
iflag=-1
eps=10.0^-24
mx=10000
xj = zeros(Float64,mx)     
fk = zeros(Complex128,mx) 

for i = 1:nj
    xj[i] = 2*pi*real(i)/ms 
end
for j = 1:round(Int,ms/2)
    fk[j] = (j+ms/2)+0.0im
end
for j = round(Int,ms/2+1):ms          
    fk[j] = (j-ms/2)+0.0im
end

cj=call_ionufft1d2(nj,ms,iflag,eps,mx,xj,fk)
println(cj)