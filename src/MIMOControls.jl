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
		complex2real

using Pkg
l = ["ControlSystems","LinearAlgebra","JuMP","SCS","Convex"]
for item in l
    Pkg.add(item)
end

using JuMP, Convex, SCS, ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
