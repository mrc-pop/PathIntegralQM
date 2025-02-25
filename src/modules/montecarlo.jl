#!/usr/bin/julia

# ------------------------------- Initialization -------------------------------

"""
Struct for a configuration with N time slicings, where SimBeta is the adimensional quantity
βε, where ε=ħ²/[m(2πR)²], and η = βε/N is the adimensional spacing.
"""
struct Configuration
    Lattice::Vector{Float64}
    SimBeta::Float64
    N::Int64
    Eta::Float64
end

"""
Create a lattice as an array of size N (number of time steps), initialized to +0.5.
"""
function SetLattice(SimBeta::Float64, N::Int64)
    Lattice = fill(0.5, N)
    Eta = SimBeta / N
    Configuration(Lattice, SimBeta, N, Eta)
end

# ------------------------------- Path processing ------------------------------ 

"""
Calculate oriented distance (diff) d(x,y), such that d(x,y) ∈ [-1/2,+1/2].
"""
function CalculateDistance(x::Float64, y::Float64)
    D = x - y
    if abs(D) <= 0.5
        return D
    elseif D > 0.5
        return D - 1.0
    else
        return D + 1.0
    end
end

"""
Precalculate neighbours. Returns two arrays: the first one contains the previous neighbours,
and the second one contains the next neighbours, both with periodic boundary conditions.
"""
function PrecalculateNeighbours(N::Int64)
    NextNeighbours = [mod(i, N) + 1 for i in 1:N]
    PrevNeighbours = [mod(i - 2, N) + 1 for i in 1:N]
    return PrevNeighbours, NextNeighbours
end

"""
Calculate winding number Q of the configuration.
"""
function CalculateQ(Config::Configuration)
    # Calculate diffs d(x_i, x_{i+1}) for all i
    Diffs = [CalculateDistance(Config.Lattice[R+1], Config.Lattice[R])
             for R in 1:(Config.N-1)]
    # Sum the diffs
    TotalDistance = sum(Diffs)
    #return round(Int, TotalDistance)
    return TotalDistance
end

# ---------------------------------- Updates -----------------------------------

"""
Single-site metropolis update. Returns the updated counter AccSteps.
"""
function MetropolisUpdate!(
    Config::Configuration,
    Site::Int64;
    Δ=0.5,
    verbose=false
)
    η = Config.Eta
    x = Config.Lattice[Site]
    xNext = Config.Lattice[Site+1]
    xPrev = Config.Lattice[Site-1]

    xTest = mod(x + Δ * (1 - 2 * rand()), 1)
	xAvg = (xNext + xPrev)/2
	
	if (abs(CalculateDistance(xAvg, xTest)) <= abs(CalculateDistance(xAvg, x)))
	 	Config.Lattice[Site] = xTest
        if verbose
        	printstyled("\nSite=$Site, x=$(round(x,digits=3)), xTest=$(round(xTest,digits=3))\n", color=:cyan)
            printstyled("ΔS<=0. Automatically accepted.\n", color=:cyan)
            
        end
        
        return 1	# Add to external counter
        
	elseif (abs(CalculateDistance(xAvg, xTest)) > abs(CalculateDistance(xAvg, x)))
		ΔS = 0.5 * η^(-1) * (
			CalculateDistance(xTest, xPrev)^2
			+ CalculateDistance(xNext, xTest)^2
			- CalculateDistance(x, xPrev)^2
			- CalculateDistance(xNext, x)^2
			)

		if verbose
			printstyled("\nSite=$Site, x=$(round(x,digits=3)), xTest=$(round(xTest,digits=3))\n", color=:yellow)
            printstyled("ΔS=$(round(ΔS,digits=3)), exp(-ΔS)=$(round(exp(-ΔS),digits=3))\n", color=:yellow)
		end

		if rand() < exp(-ΔS)
			Config.Lattice[Site] = xTest
			if verbose
				printstyled("Accepted. Yay.\n", color=:green)
			end
			
			return 1	# Add to external counter
		else
			if verbose
			    printstyled("\nNot accepted. Rip.\n", color=:red)
			end
			
			return 0	# Add to external counter
		end
	end
end

"""
Heatbath update.
"""
function HeatBathUpdate!(
	Config::Configuration,
	Site::Int64;
	verbose=false
)
    # Gaussian parameter linked to simulations parameters
    Alpha = 1/(Config.Eta * Config.SimBeta)

    # Interpolate lowest action postion
    Center = (Config.Lattice[Site-1] + Config.Lattice[Site+1])/2

    # Extract update from gaussian pdf with null center e unitary variance,
    # then recenter the result
    Update = mod(Center + randn()/(2*sqrt(Alpha)), 1)
    Config.Lattice[Site] = Update
    
    return
end

"""
Cluster \"tailor\" update [Bonati, D'Elia 2018], with tolerance ε.
"""
function TailorUpdate!(Config::Configuration, Site::Int64, ε::Float64)
    i = Site
    x = Config.Lattice[i]
    iEnd = findfirst(y -> abs.(CalculateDistance(x + 0.5, y))<ε, Config.Lattice)

    # Propose path where y ↦ x_0 - (y - x_0) = 2 x_0 - y
    # (inversion about x_0) for y in {x_i+1, ..., x_iEnd}
    xxTest = [mod(2 * x - Config.Lattice[j], 1) for j in (i+1):iEnd]

    xxNewPlot = vcat(Config.Lattice[i], xxTest, Config.Lattice[iEnd+1])

    # Calculate change in action after update
    ΔS = 0.5 * Config.Eta^(-1) * (
        CalculateDistance(Config.Lattice[iEnd+1], xxTest[end])^2
        -CalculateDistance(Config.Lattice[iEnd+1], Config.Lattice[iEnd])^2
        )

    if ΔS <= 0
        Config.Lattice[(i+1):iEnd] = xxTest # accept the update
    elseif rand() < exp(-ΔS)
        Config.Lattice[(i+1):iEnd] = xxTest # accept the update
    end

    return iEnd, xxNewPlot
end

# ------------------------------------ Main ------------------------------------

function main()
    NMetro = 100000
    N = 20

    Config = SetLattice(1.0, N)
    Counter = 0

    println("Performing $NMetro Metropolis steps...")

    for i in 1:NMetro
        # Site = i % (N-2) + 2 # from 2 to N-1 (sequential)
        Site = rand(2:N-1) # (random)
        Counter += MetropolisUpdate!(Config, Site; Δ=0.5, verbose=true)
        # HeatBathUpdate!(Config, Site)
    end

    println("Accepted steps: $Counter/$NMetro")

    p = PlotPath(Config; location="nw")

    # println("Performing tailor update...")
    # i = 2
    # ε = 0.1
    # iEnd, xxNewPlot = TailorUpdate!(Config, i, ε)

    # PlotTailorUpdate!(Config, i, iEnd, xxNewPlot)

end

if abspath(PROGRAM_FILE) == @__FILE__

	PROJECT_ROOT = @__DIR__ # Absloute path up to .../PathIntegralQM/src

	include(PROJECT_ROOT * "/../setup/graphic_setup.jl")
	include(PROJECT_ROOT * "/plots.jl")

	main()
end
