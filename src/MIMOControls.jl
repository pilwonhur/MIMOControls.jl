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
pkg.add(JuMP)
pkg.add(SCS)
pkg.add(Convex)
pkg.add(ControlSystems)
pkg.add(LinearAlgebra)

using JuMP, Convex, SCS, ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
