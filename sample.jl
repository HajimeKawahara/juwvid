#!/Applications/Julia-0.4.2.app/Contents/Resources/julia/bin/julia 
using ArgParse
using Winston
using cohenclass

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-m"
            help = "Choose the method (wv=Wigner Ville, ...)"
            required = true
    end
    return parse_args(s)
end

function main()
    println(PyDict(pyimport("matplotlib")["rcParams"])["backend"])
    parsed_args = parse_commandline()
    method=println(parsed_args["m"])
    
    nsamp=100
    x=linspace(0.0,100pi,nsamp)
    wp=1.0 #angular frequency of propagation
    wf=0.01 #modulation
    s=1.0 # sterngth
    y=cos(wp*x+s*sin(wf*x))
    display(plot(x, y))

end

main()