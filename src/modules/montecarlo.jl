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
    N = Config.N
    xNext = Config.Lattice[mod1(Site+1,N)]
    xPrev = Config.Lattice[mod1(Site-1,N)] # 1 ↦ N

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

	N = Config.N	
	
    # Gaussian parameter linked to simulations parameters
    Alpha = 1/(Config.Eta * Config.SimBeta)

    # Interpolate lowest action postion
    Center = (Config.Lattice[mod1(Site-1,N)] + Config.Lattice[mod1(Site+1,N)])/2

    # Extract update from gaussian pdf with null center e unitary variance,
    # then recenter the result
    Update = mod(Center + randn()/(2*sqrt(Alpha)), 1)
    Config.Lattice[Site] = Update

    return
end

"""
Cluster \"tailor\" update [Bonati, D'Elia 2018], with tolerance ε.
"""
function TailorUpdate!(Config::Configuration, Site::Int64, ε::Float64; verbose=false)

    i = Site
    x = Config.Lattice[i]

    # Shift lattice, so to have Lattice[i+1] in first position
    circshift!(Config.Lattice, -i)

    # Find iEnd as prescribed by the tailor algorithm
    DeltaiEnd = findfirst(y -> abs.(CalculateDistance(x + 0.5, y))<ε, Config.Lattice)

    if DeltaiEnd === nothing
        verbose && println("No iEnd found")
        return [0, 0, nothing, nothing]
    else
        Found = 1
    end

    # Propose path where y ↦ x_0 - (y - x_0) = 2 x_0 - y
    # (inversion about x_0) for y in {x_i+1, ..., x_iEnd}
    xxTest = [mod(2 * x - Config.Lattice[j], 1) for j in 1:DeltaiEnd]

    # Calculate change in action after update
    ΔS = 0.5 * Config.Eta^(-1) * (
        CalculateDistance(Config.Lattice[mod1(DeltaiEnd+1,N)], xxTest[end])^2
        -CalculateDistance(Config.Lattice[mod1(DeltaiEnd+1,N)], Config.Lattice[DeltaiEnd])^2
        )

    Acc = 0

    if ΔS <= 0 || rand() < exp(-ΔS)
        Config.Lattice[1:DeltaiEnd] = xxTest # accept the update
        Acc = 1
        verbose && println("Tailor accepted")
    else
        verbose && println("Tailor rejected")
    end

    # Shift back lattice
    circshift!(Config.Lattice, i)
    iEnd = mod1(i + DeltaiEnd,N)

    xxNewPlot = vcat(Config.Lattice[i], xxTest, Config.Lattice[mod1(iEnd+1,N)])

    return Found, Acc, iEnd, xxNewPlot
end

# ------------------------------------ Main ------------------------------------

const NSweepsTherm = 10000
const NSweeps = Int(1e6)               # number of updates of the whole lattice
const Δ = 0.3
const sequential = false
const NN = [100, 200, 300]                    # number of lattice points
const SimBeta = 2.0
const εε_over_ηη = [0.5, 1.0, 1.5, 2.0]      # tolerance of tailor method (article)

using Statistics

function main()

    for N in NN
        printstyled("\nWorking on N=$N", color=:yellow)
        η = SimBeta/N

        for Fraction in εε_over_ηη
            ε = Fraction * η
            println("\nUsing ε=$ε, N=$N, η=$η, ε/η = $Fraction")

            Config = SetLattice(1.0, N)
            FoundTailor = 0
            AccTailor = 0

            for _ in 1:(NSweepsTherm*N)
                Site = rand(2:N-1)
                MetropolisUpdate!(Config, Site; Δ)
            end

            LL = fill(0, NSweeps)

            for i in 1:NSweeps
                # Sweep over lattice
                for j in 1:N
                    Site = rand(2:N-1)
                    MetropolisUpdate!(Config, Site; Δ)
                end

                Site = rand(2:Int(N/2))
                Found, Acc, iEnd, _ = TailorUpdate!(Config, Site, ε; verbose=false)
                if iEnd !== nothing
                    LL[i] = iEnd - Site
                end

                FoundTailor += Found
                AccTailor += Acc
            end

            # exclude zeros
            LL = LL[LL .!= 0]

            MeanL = mean(LL)
            StdL = std(LL)

            println("Over $NSweeps tailor updates, I found iEnd $FoundTailor times (ratio $(round(FoundTailor/NSweeps,digits=4)))")
            println("and accepted the tailor metropolis step $AccTailor times (acceptance $(round(AccTailor/FoundTailor, digits=4)))")
            println("Length of the cluster: L=$(round(MeanL, digits=2)) ± $(round(StdL,digits=2))")
        end
    end
end

    # println("Performing $NMetro Metropolis steps...")

    # for i in 1:NMetro
    #     # Site = i % (N-2) + 2 # from 2 to N-1 (sequential)
    #     Site = rand(2:N-1) # (random)
    #     Counter += MetropolisUpdate!(Config, Site; Δ=0.5, verbose=true)
    #     # HeatBathUpdate!(Config, Site)
    # end

    # println("Accepted steps: $Counter/$NMetro")

    # p = PlotPath(Config; location="nw")

    # # println("Performing tailor update...")
    # # i = 2
    # # ε = 0.1
    # # iEnd, xxNewPlot = TailorUpdate!(Config, i, ε)

    # # PlotTailorUpdate!(Config, i, iEnd, xxNewPlot)

if abspath(PROGRAM_FILE) == @__FILE__

	PROJECT_ROOT = @__DIR__ # Absloute path up to .../PathIntegralQM/src

	include(PROJECT_ROOT * "/plots.jl")

	main()
end
