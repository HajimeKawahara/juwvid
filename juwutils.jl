module juwutils

function index_to_frequency(indf, fin, dx, nsample)
    # indf : indices
    # fin : 
    # dx : the size of the time bin (x-axis)
    # nsample : # of data
    if isnan(fin[1])
        fin=collect(1:length(indf))
    end

    offset=(fin[end]-fin[1])/nsample
    freqfac=1/nsample/dx/2
    return (fin[round(Int,indf)]-offset)*freqfac    
end

end