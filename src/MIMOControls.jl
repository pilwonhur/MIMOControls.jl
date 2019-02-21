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
		hinflmi,
		h2lmi,
		splitSS



using JuMP, SCS, ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
