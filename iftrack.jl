module iftrack
# Instantaneous Frequency Tracking

function iterate_mcif(timeseq,htwindow,hfwindow,iguess,tfr,prec)
    iseq=Int64[]
    nf=size(tfr)[1]
    nt=size(tfr)[2]
    for itime=timeseq
        ti=max(1,itime-htwindow)
        te=min(nt,itime+htwindow)
        ct=converge_center2d(iguess,hfwindow,abs(vec(tfr[:,ti:te])),nf,prec)
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
        fe=min(center_test+hsize,nend)
        fs=max(1,center_test-hsize)
        pos=collect(fs:fe)
        tip=abs(array2d[fs:fe,:])
        ceni=round(Int,sum(pos.*tip)/sum(tip))
        diff=abs(center_test-ceni)
        #print(nstep,":",diff,"\n")
        center_test=ceni
        nstep=nstep+1
    end
    return center_test
end

end
