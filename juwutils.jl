module juwutils

function index_to_frequency(indf, fin, dx)
    # indf : indices
    # fin : 
    # dx : the size of the time bin (x-axis)
    if isnan(fin[1])
        fin=collect(1:length(indf))
    end

    nsample=length(fin)
    offset=(fin[end]-fin[1])/nsample
    freqfac=1/nsample/dx/2
    return (fin[round(Int,indf)]-offset)*freqfac    
end

end