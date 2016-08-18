module juwutils

function index_to_frequency(indf, fin, dx, nsample, nft=NaN, finend=NaN,fin1=NaN, Nz=NaN)
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
        if isnan(fin1) fin1 = 1.0 end
        if isnan(finend) finend= nsample end
        if isnan(Nz)
            offset=(finend-fin1)/nsample
            return (indf-offset)*freqfac    
        else 
            offset=(finend-fin1)/nsample
            return (indf+round(Int,fin1*Nz)-1-offset)*freqfac/Nz

        end
    else
#        offset=(fin[end]-fin[1])/nsample
        offset=(fin[end]-fin[1])/nft
        return (fin[round(Int,indf)]-offset)*freqfac    
    end
end

function frequency_to_index(freqarray, dx, nsample, nft=NaN, df=NaN)
    # dx : the size of the time bin (x-axis)
    # nsample : # of data
    # nft : # of frequency bins
    if isnan(nft)
        nft=nsample
        println("Assuming nft = nsample.")
    end

#    if isnan(df) && size(freqarray)>2
#        println("Warning: Assuming |f[2] - f[1]| as df.")
#        offset=abs(freqarray[2]-freqarray[1])/nft
#    elseif isnan(df)
#        println("Error: Provide df.")
#    else
#        offset=df/nft
#    end
    offset=0.0
    freqfac=1/nsample/dx/2
    return freqarray/freqfac + offset + 1
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

function jbinning(x,y,medbinsize=12,mode="median")
yb=Float64[]
xb=Float64[]
for i=1:round(Int,(length(y)-1)/medbinsize)-1    
    k=(i-1)*medbinsize+1
    l=i*medbinsize
    if mode=="median"
        xbtmp=median(x[k:l])
#        xtmp=median(itc[k:l]*dt+start)
        ybtmp=median(y[k:l])
    elseif mode=="mean"
        xbtmp=mean(x[k:l])
        ybtmp=mean(y[k:l])
    end
    append!(xb,[xbtmp])
    append!(yb,[ybtmp])
end

return xb, yb
end

end
