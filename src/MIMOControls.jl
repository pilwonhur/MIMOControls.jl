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
		daug



using ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
