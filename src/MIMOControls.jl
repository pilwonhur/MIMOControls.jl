module MIMOControls

export	kalmandecomp,
		findprojector,
		mqr,
		mmult,
		madd,
		msub,
		transp,
		cjt,
		cj,
		minv,
		mscl,
		sbs,
		abv,
		daug,
		minimumreal,
		eye,
		eyec,
		hinflmi,
		hinfbis,
		h2lmi,
		h2gram,
		splitSS,
		complex2real,
		LQRlmi,
		h2syn,
		hinfsyn,
		lft,
		are

using Pkg
l = ["ControlSystems","LinearAlgebra","JuMP","ProxSDP"]
for item in l
    Pkg.add(item)
end

using JuMP, ProxSDP, ControlSystems, LinearAlgebra
# using SCS, Convex

include("matrix_comps.jl")
end # module
