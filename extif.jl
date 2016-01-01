module extif
#extract instant frequency

function maxif(tfr)
    nf,nt=size(tfr)
    indf=zeros(nt)
    for it=1:nt
        i=indmax(tfr[:,it])
        indf[it]=real(i)
    end 
    
    return indf
end

function aveif(tfr,isf=NaN,ief=NaN)
    nf,nt=size(tfr)
    if isnan(isf) isf=1 end 
    if isnan(ief) ief=nf end 
    
    nnf=ief-isf+1
    indf=zeros(nt)
    ifreq=collect(isf:ief)
    for it=1:nt
        i=sum(ifreq[:].*tfr[isf:ief,it])/sum(tfr[isf:ief,it])
        indf[it]=i
    end 
    
    return indf
end


end