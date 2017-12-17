module iftrack
import juwutils

# Instantaneous Frequency Tracking

function iterate_mcif(timeseq,htwindow,hfwindow,iguess,tfr,prec=1)
    iseq=Int64[]
    nf=size(tfr)[1]
    nt=size(tfr)[2]
    for itime=timeseq
        ti=max.(1,itime-htwindow)
        te=min.(nt,itime+htwindow)
        ct=converge_center2d(iguess,hfwindow,abs.(vec(tfr[:,ti:te])),nf,prec)
        iguess=ct
        push!(iseq,iguess)
        #print(itime,"-",ct,"\n")
    end
    return collect(timeseq), iseq
end

function converge_center2d(inguess,hsize,array2d,nend,prec=1,nlim=1000)
    #hsize: half size of frequency window
    #nlim: iteration limit
    n=2*hsize+1
    center_test=inguess
    diff=prec+1
    nstep=1
    while diff>prec && nstep < nlim
        fe=min.(center_test+hsize,nend)[1]
        fs=max.(1,center_test-hsize)[1]
        pos=collect(fs:fe)
        tip=abs.(array2d[fs:fe,:])
        ceni=round.(Int,sum(pos.*tip)/sum(tip))
        diff=abs.(center_test-ceni)[1]
        #print(nstep,":",diff,"\n")
        center_test=ceni
        nstep=nstep+1
    end
    return center_test
end

function track_mode(tfr,t,tinit,finit,tend,windowdf,tfrmode=2,rev=0,nthin=10,htwindow=1)
    #tfrmode=2: STFT,SM,LW,PolyWV
    nf=size(tfr)[1]
    nsample=size(tfr)[2]
    dt=t[2]-t[1]
    itstart=max(round(Int,(tinit-t[1])/dt),1)
    itend=min(round(Int,(tend-t[1])/dt),nsample)
    
    iguess=round.(Int,juwutils.frequency_to_index([finit], dt, nsample)/tfrmode)[1]
    wguess=round.(Int,juwutils.frequency_to_index([finit+windowdf], dt, nsample)/tfrmode)[1]
    hfwindow=round(Int,wguess-iguess)
    timeseq=itstart:nthin:itend #sparse sampling
    if rev==1
        timeseq=reverse(timeseq)
    end
    timeseq,iseq=iftrack.iterate_mcif(timeseq,htwindow,hfwindow,iguess,tfr);
    fx=juwutils.index_to_frequency(iseq.*tfrmode,NaN,dt,nsample);
    return t[timeseq],fx
end

function get_poltrack(timeseq,iseq,I,Q,U,V)
    #polarization tracking
    #tip=3
    nv=size(I)[1]
    n=length(iseq)
    Ival=zeros(n);Qval=zeros(n);Uval=zeros(n);Vval=zeros(n);
    for i=1:n
        it=timeseq[i]
        jt=iseq[i]
        Ival[i]=I[jt,it]
        Qval[i]=Q[jt,it]
        Uval[i]=U[jt,it]
        Vval[i]=V[jt,it]
        #for kk=maximum([jt-tip,1]):minimum([nv,jt+tip])
        #    Ival[i]+=I[kk,it]
        #    Qval[i]+=Q[kk,it]
        #    Uval[i]+=U[kk,it]
        #    Vval[i]+=V[kk,it]
        #end
    end
    return Ival,Qval,Uval,Vval
end

end

