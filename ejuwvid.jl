println("Starting interactive juwvid...")
println("Julia Version=",VERSION)
if VERSION < v"0.6-"
    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    println("You use unsupported version of julia for juwvid. Update julia.")
    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
end
file=ARGS[1]
#method=ARGS[1]
data=readdlm(file)
t=data[:,1]
x=data[:,2]

using PyPlot
PyPlot.plot(t,x)
PyPlot.xlabel("t")
PyPlot.ylabel("data")
PyPlot.show()

import cohenclass
import estif
import extif
import jnufft
import juwutils
import lwigner
import polywv
import pwvaw
import smethod
import stft
import juwplot
import pm
import iftrack

