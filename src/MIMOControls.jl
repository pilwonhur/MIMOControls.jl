module MIMOControls

export	kalmandecomp,
		findprojector

using ControlSystems, LinearAlgebra

include("matrix_comps.jl")
end # module
