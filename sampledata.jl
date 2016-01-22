module sampledata

function genfm(nsamp,wp=1.0,wf=0.01,s=1.0, xend=100*pi)
    x=linspace(0.0,xend,nsamp)
    y=cos(wp*x+s*sin(wf*x))
    iw=wp+s*wf*cos(wf*x)
    ynorm=pi/x[end]
    return x, y, iw, ynorm
end

function genlinfm(nsamp,wp=1.0,s=0.01, xend=100*pi)
    x=linspace(0.0,xend,nsamp)
    y=cos(wp*x+0.5*s*x.*x)
    iw=wp+s*x
    ynorm=pi/x[end]
    return x, y, iw, ynorm
end

function genstepfm(nsamp,xend=1.0)
    nsamp=1024
    x=linspace(0.0,xend,nsamp)
    iw=nsamp/4*atan(250*(x-0.5))+nsamp/4
    phase=cumsum(iw)/nsamp
    y=cos(phase)
    ynorm=pi/x[end]
    return x, y, iw, ynorm
end


end