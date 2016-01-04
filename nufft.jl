#### DO NOT WORK YET !!!!!! ####

module nufft

#type 2
function getl(i,N)
    l=round(Int64,rem(N+i,N)+1)
    return l
end

function next235(base,limit=10^12)
    #limit : maximum value for search
    if rem(base,2)==0
        nextv=base
    else
        nextv=base+1
    end
    p=0
    q=0
    r=0
    while true
        val=nextv
        while true
            if rem(round(Int,val),2) != 0 break end
            p +=1 
            val=round(Int,val/2)  
        end
        #println(val," p:",p)
        while true
            if rem(round(Int,val),3) != 0 break end
            q +=1 
            val=round(Int,val/3)  
        end    
        #println(val," q:",q)
        while true
            if rem(round(Int,val),5) != 0 break end
            r +=1 
            val=round(Int,val/5)  
        end    
        #println(val," r:",r)

        if val == 1 break end
        nextv += 2
        p=0
        q=0
        r=0
        if nextv > limit 
            println("Exceed Maximum value=", limit)
            nextv=0
            break 
        end 
    end
    return nextv, p, q, r
end

function nufft1d2(nj,M,xj,fk,iflag=-1,eps=1.e-24)
    #fw 0:iwtot -> + 1 real
    #fk use getl
    #xc use getl
    nxc=147*2+1
    xc=zeros(nxc)

    if eps <= 1E-33 || eps >= 1E-1
        println("eps =",eps,"must satisfy 1e-33 < eps < 1e-1.")
    end
    if eps > 1E-11 
        ratio = 2 
    else 
        ratio = 3
    end

    nspread = round(Int64, -log(eps)/(pi*(ratio - 1)/(ratio - 0.5)) + 0.5)
    nf1 = ratio*M
    if 2*nspread > nf1 nf1,p1,q1,r1 = next235(2*nspread) end 

    rrlamb = ratio*ratio*nspread/(ratio*(ratio - 0.5))
    hx = 2*pi/nf1


#     ---------------------------------------------------------------
#     Precompute spreading constants and initialize fw
#     to hold one term needed for fast Gaussian gridding 
#     ---------------------------------------------------------------
    iw1 = 2*nf1
    iwsav = iw1 + nspread + 1
    iwtot = iwsav + 4*nf1 + 15    
#    allocate ( fw(0:iwtot))
    fw=zeros(Real,iwtot+1) #fw: i->i+1

    t1 = pi/rrlamb
    for k1 = 1:nspread
        fw[iw1+k1+1] = exp(-t1*k1^2)
    end
#    call zffti(nf1,fw(iwsav))
#
#     ---------------------------------------------------------------
#     Deconvolve and compute inverse 1D FFT
#     (A factor of (-1)**k is needed to shift phase.)
#     ---------------------------------------------------------------
#
    t1 = pi*rrlamb/nf1^2
    cross1 = 1.0/sqrt(rrlamb)
    zz = cross1*fk[1]
    fw[1] = real(zz)
    fw[2] = imag(zz)
    for k1=1:(M-1)/2
        cross1 = -cross1
        cross = cross1*exp(t1*k1^2)        
        zz = cross*fk[round(Int,k1+1)]
        fw[round(Int,2*k1+1)] = real(zz)
        fw[round(Int,2*k1+2)] = imag(zz)
        l=getl(-k1,M)
        println(l)
        zz = cross*fk[l]
        l=getl(2*(nf1-k1),M)
        fw[l] = real(zz)
        l=getl(2*(nf1-k1)+1,M)
        fw[l] = imag(zz)
    end
    cross = -cross1*exp(t1*(M/2)^2)
    if M == nextpow2(M)
        l=getl(-M/2,M)
	zz = cross*fk[l]
        l=getl(2*nf1-M,M)
        fw[l] = real(zz)
        l=getl(2*nf1-M+1,M)
        fw[l] = imag(zz)
    end
    for k1=(M+1)/2:nf1-M/2-1
        fw[round(Int,2*k1+1)] = 0.0
        fw[round(Int,2*k1+2)] = 0.0
    end
    
      #complex fft
    ffw=fft(fw[1:2*M])

    fw[1]=0.0
    fw[2]=0.0
    for k1=1:M-1        
        fw[round(Int,2*k1+1)] = real(ffw[k1])
        fw[round(Int,2*k1+2)] = imag(ffw[k1])
    end

      #call zfftf(nf1,fw(0),fw(iwsav))

#     ---------------------------------------------------------------
#     Loop over target points (1,...,nj)
#
#       1. find closest mesh point (with periodic wrapping if needed)
#       2. get contributions from regular fine grid to target
#          locations using Gaussian convolution.
#     ---------------------------------------------------------------
    t1 = pi/rrlamb
    cj=zeros(Complex64,nj)

    for j = 1:nj
        cj[j] = 0.0+0.0im
        jb1 = round(Int,(xj[j]+pi)/hx)
        diff1 = (xj[j]+pi)/hx - jb1
#       
        jb1 = rem(jb1, nf1)
        if jb1 < 0; jb1=jb1+nf1 end
        xc[1] = exp(-t1*diff1^2)
        cross = xc[1]
        cross1=exp(2*t1*diff1)
       for k1 = 1:nspread
            cross = cross * cross1
            xc[k1+1] = fw[iw1+k1+1]*cross
        end
        cross = xc[1]
        cross1 = 1.0/cross1

        for k1 = 1:nspread-1
            cross = cross * cross1
            l=getl(-k1,nxc)           
            xc[l] = fw[iw1+k1+1]*cross
        end

        jb1d = minimum([nspread-1, jb1])
        jb1u = minimum([nspread, nf1-jb1-1])
        for k1 = -nspread+1:-jb1d-1
	    zz = fw[2*(jb1+k1+nf1)+1]+fw[2*(jb1+k1+nf1)+2]*im
            l=getl(k1,nxc)
            cj[j] = cj[j] + xc[l]*zz
        end
#
        for k1 = -jb1d:jb1u
	    zz = fw[2*(jb1+k1)+1]+fw[2*(jb1+k1)+2]*im
            l=getl(k1,nxc)
            cj[j] = cj[j] + xc[l]*zz ####
        end
#

        for k1 = jb1u+1:nspread
	    zz = fw[2*(jb1+k1-nf1)+1]+fw[2*(jb1+k1-nf1)+2]*im
            l=getl(k1,nxc)
            cj[j] = cj[j] + xc[l]*zz
        end
        println("-")
 

    end
    
    return cj
end

end

#import DSP
M=8
nj=8
i=1
xj=zeros(nj)
fk0=zeros(Complex64,M)
for i = 1:nj 
    xj[i] = 2*pi*real(i)/M
end
for j = 1:round(Int,M/2)
    fk0[j]= real(j+M/2)+0.0im
end
for j = round(Int,M/2+1):M
    fk0[j]= real(j-M/2)+0.0im
end
cj=nufft.nufft1d2(nj,M,xj,fk0)
println(cj)

########################################################
# Note 
#
# nufft.jl was imported from NUFFT Software with modification.
#
# information on NUFFT: 
#Copyright (c) 2009-2013, Leslie Greengard, June-Yub Lee and Zydrunas Gimbutas
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met: 

#1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer. 
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution. 

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#The views and conclusions contained in the software and documentation are those
#of the authors and should not be interpreted as representing official policies, 
#either expressed or implied, of the FreeBSD Project.

