#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/modules/montecarlo.jl")

# Import MetropolisUpdate!, HeatBathUpdate!, TailorUpdate!

function main()
	print("Hello World!")
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
