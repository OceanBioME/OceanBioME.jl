using Profiling # Prepend something to distinguish it from build-in Profile package
using Profile
# using PProf
# using ProfileCanvas
using ProfileSVG

@info "Precompiling"
@time precompile(Profiling.Cases.simple_LOBSTER, ())
@info "END"

@info "Fast Kill run - Does more precompilation"
@time Profiling.Cases.simple_LOBSTER(;fast_kill=true)
@info "END"


Profile.clear()

Profile.init(n = 10^8, delay=0.001)
@info "True run"
@time @profile Profiling.Cases.simple_LOBSTER()
@info "END"


open("prof.txt","w") do io
    Profile.print(io, noisefloor=2.0)
end

ProfileSVG.save( "prof.svg")
