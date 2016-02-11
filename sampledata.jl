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

function genstepfm(nsamp=1024,xend=1.0)
    x=linspace(0.0,xend,nsamp)
    iw=nsamp/4*atan(250*(x-0.5))+pi*nsamp/4
    phase=cumsum(iw)/nsamp
    y=cos(phase)
    iy=sin(phase)*im
    ynorm=pi/x[end]
    return x, y+iy, iw, ynorm
end

function genmultifm622(nsamp=512)
    #Boashash+2015 Example 6.2.2
    t=collect(linspace(-1.0,1.0,nsamp))
    x=exp(-t.*t).*cos(25.0*pi*t)+cos(120.0*t.*t.*t+45.0*pi*t)+1.5*exp(-25.0*t.*t).*cos(40.0*pi*t.*t+150.0*pi*t)
    return t, x
end

function genmultifm622x(nsamp=512)
    #Boashash+2015 modified Example 6.2.2
    t=collect(linspace(-1.0,1.0,nsamp))
    x=exp(-t.*t).*cos(25.0*pi*t)+cos(12.0*t.*t.*t+40.0*pi*t)+1.5*exp(-25.0*t.*t).*cos(4.0*pi*t.*t+65.0*pi*t)
    return t, x
end

function genmultifm623(nsamp=256)
    #Boashash+2015 Example 6.2.3
    t=collect(linspace(-1.0,1.0,nsamp))
    x=cos(20.0*sin(pi*t)+30.0*pi*t)+sin(20.0*cos(pi*t)+100.0*pi*t)
    return t, x
end

function genmultifm623x(nsamp=256)
    #Boashash+2015 Example 6.2.3
    t=collect(linspace(-1.0,1.0,nsamp))
    x=cos(2.0*sin(pi*t)+30.0*pi*t)+sin(2.0*cos(pi*t)+60.0*pi*t)
    return t, x
end

end

