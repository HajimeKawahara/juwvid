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

end