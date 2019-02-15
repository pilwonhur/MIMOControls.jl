module MIMOControls

export	greet,
		greet2

greet() = print("Hello World!")
greet2() = print("Hello World 2 !")

using ControlSystems, Polynomials, DifferentialEquations, LinearAlgebra

include("matrix_comps.jl")
end # module
