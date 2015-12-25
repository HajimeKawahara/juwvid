module sampledata

function genfm(nsamp,wp=1.0,wf=0.01,s=1.0)
    x=linspace(0.0,100pi,nsamp)
    y=cos(wp*x+s*sin(wf*x))
    return x, y
end

end