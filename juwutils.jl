module juwutils

function index_to_frequency(indf, fin, dx, nsample, nft=NaN, finend=NaN,fin1=NaN)
    # indf : indices
    # fin : 
    # dx : the size of the time bin (x-axis)
    # nsample : # of data
    # nft : # of frequency bins

#    if isnan(fin[1])
#        fin=collect(1:length(indf))
#    end
    if isnan(nft)
        nft=nsample
        println("Assuming nft = nsample.")
    end

    freqfac=1/nsample/dx/2
    if isnan(fin[1])
        offset=(finend-fin1)/nsample
        return (indf-offset)*freqfac    
    else
#        offset=(fin[end]-fin[1])/nsample
        offset=(fin[end]-fin[1])/nft
        return (fin[round(Int,indf)]-offset)*freqfac    
    end
end

function frequency_to_index(freqarray, dx, nsample, nft=NaN)
    # dx : the size of the time bin (x-axis)
    # nsample : # of data
    # nft : # of frequency bins

    if isnan(nft)
        nft=nsample
        println("Assuming nft = nsample.")
    end

    offset=abs(freqarray[2]-freqarray[1])/nft
    freqfac=1/nsample/dx/2
#    return round(Int,freqarray/freqfac + offset)
    return freqarray/freqfac + offset
end


function degradexy(x,y,smoothf) 
    nsample=length(y)
    nT=round(Int,nsample/smoothf)
    meanx=zeros(nT)
    meany=zeros(nT)
    for i=1:nT
        meanx[i]=mean(x[(i-1)*smoothf+1:i*smoothf])
        meany[i]=mean(y[(i-1)*smoothf+1:i*smoothf])
    end
    return meanx,meany
end

function throw_intarray(t,lc,dt,offt,fillval=1.0)
    #dt=t[2]-t[1]
    #offt=t[1]
    jend=round(Int,(t[end]-offt)/dt+1.0)
    lcn=ones(jend)*fillval;
    for i=1:length(t)
        if !isnan(t[i])    
            x=(t[i]-offt)/dt+1.0
            j=round(Int,x)
            lcn[j]=lc[i]            
    #println(i,"-",j)
        end
    end
    return lcn
end

end