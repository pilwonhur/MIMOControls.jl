module MIMOControls

export	kalmandecomp,
		findprojector,
		mqr,
		mmult,
		madd,
		msub,
		transp,
		cjt,
		minv,
		mscl,
		sbs,
		abv,
		daug,
		minimumreal,
		eye



using JuMP, SCS, ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
