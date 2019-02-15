module MIMOControls

export	kalmandecomp,
		findprojector,
		mqr

using ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
