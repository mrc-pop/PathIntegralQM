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
function CalculateQ(Config::Configuration; WarningQty=0.001)
    # Calculate diffs d(x_i, x_{i+1}) for all i
    Diffs = [CalculateDistance(Config.Lattice[mod1(R+1,Config.N)], Config.Lattice[R])
             for R in 1:Config.N]
    # Sum the diffs
    TotalDistance = sum(Diffs)

    if abs(TotalDistance - round(TotalDistance)) > WarningQty
        @warn "Q is more than $WarningQty away from an integer"
        return TotalDistance
    end

    return round(Int64, TotalDistance)
end

"""
Calculate Euclidean action S of the configuration.
"""
function CalculateS(Config::Configuration)
    # Calculate diffs d(x_i, x_{i+1}) for all i
    Diffs = [CalculateDistance(Config.Lattice[mod1(R+1,Config.N)], Config.Lattice[R])
             for R in 1:Config.N]
    # Sum the diffs
    S = 0.5 * (Config.Eta)^(-1) * sum(Diffs.^2)

    return S
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
    xNext = Config.Lattice[mod1(Site+1,N)] # N ↦ 1
    xPrev = Config.Lattice[mod1(Site-1,N)] # 1 ↦ N

    xTest = mod(x + Δ * (1 - 2 * rand()), 1)

    ΔS = 0.5 * η^(-1) * (
        CalculateDistance(xTest, xPrev)^2
        + CalculateDistance(xNext, xTest)^2
        - CalculateDistance(x, xPrev)^2
        - CalculateDistance(xNext, x)^2
        )

    if verbose
        printstyled("\nSite=$Site, x=$(round(x,digits=3)), xTest=$(round(xTest,digits=3))\n",
            color=:yellow)
        printstyled("ΔS=$(round(ΔS,digits=3)), exp(-ΔS)=$(round(exp(-ΔS),digits=3))\n",
            color=:yellow)
    end

    if ΔS <= 0 || rand() < exp(-ΔS)
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

"""
Heatbath update.
"""
function HeatBathUpdate!(
	Config::Configuration,
	Site::Int64;
	verbose=false
)

	N = Config.N

    xNext = Config.Lattice[mod1(Site+1,N)] # N ↦ 1
    xPrev = Config.Lattice[mod1(Site-1,N)] # 1 ↦ N

    # Gaussian parameter linked to simulations parameters
    Alpha = 1/(Config.Eta * Config.SimBeta)
    # OLD: # Center = (Config.Lattice[mod1(Site-1,N)] + Config.Lattice[mod1(Site+1,N)])/2

    # Interpolate lowest action postion
    Dist = CalculateDistance(xNext, xPrev)
    Center = mod(xPrev + Dist/2, 1)

    # Extract update from gaussian pdf with null center e unitary variance,
    # then recenter the result
    Update = mod(Center + randn() * sqrt(Config.Eta/2), 1)
    Config.Lattice[Site] = Update

    return
end

"""
Cluster «tailor» update [Bonati, D'Elia 2018], with tolerance ε.
"""
function TailorUpdate!(Config::Configuration, Site::Int64, ε::Float64; verbose=false)

    # Simplify notation
    i = Site
    x = Config.Lattice[i]
    N = Config.N

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

    if verbose
        printstyled("\nSite=$Site, x_S=$(round(x,digits=3)), x_E=$(round(Config.Lattice[iEnd],digits=3))\n",
            color=:yellow)
        printstyled("ΔS=$(round(ΔS,digits=3)), exp(-ΔS)=$(round(exp(-ΔS),digits=3))\n",
            color=:yellow)
    end

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


"""
Perform a parallel tempering update, proposing to swap Config1 and Config2.
"""
function ParallelTemperingUpdate!(
    Config1::Configuration,
    Config2::Configuration;
    verbose=false
)

    if length(Config1.Lattice) != length(Config2.Lattice)
        error("Error! Lattices must be of the same size to propose swap.")
    end

    S1 = CalculateS(Config1)
    η1 = Config1.Eta

    S2 = CalculateS(Config2)
    η2 = Config2.Eta

    # Calculate the new action after the exchange of S1 ↔ S_2
    S1New = η2 / η1 * S2
    S2New = η1 / η2 * S1

    # Define acceptance probability
    AccProb = exp(- S2New - S1New + S1 + S2)

    if verbose
        println("Acceptance probability of swap: $AccProb.")
    end

    # Metropolis acceptance
    if AccProb > 1 || rand() < AccProb
        # Accept swap (using a temporary array)
        temp = copy(Config1.Lattice)
        Config1.Lattice .= Config2.Lattice
        Config2.Lattice .= temp
        if verbose
            println("Swap accepted.")
        end
        return 1
    end

    if verbose
        println("Swap rejected.")
    end

    return 0
end
