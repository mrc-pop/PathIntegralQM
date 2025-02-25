#!/usr/bin/julia

include("graphic_setup.jl")

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

"""
Single-site metropolis update. Returns the updated counter AccSteps.
"""
function MetropolisUpdate!(
    Config::Configuration,
    Site::Int64,
    AccSteps::Int64;
    Δ=0.5,
    verbose=false
)
    η = Config.Eta
    x = Config.Lattice[Site]
    xNext = Config.Lattice[Site+1]
    xPrev = Config.Lattice[Site-1]

    xTest = mod(x + Δ * (1 - 2 * rand()), 1)

    ΔS = 0.5 * η^(-1) * (
        CalculateDistance(xTest, xPrev)^2
        + CalculateDistance(xNext, xTest)^2
        - CalculateDistance(x, xPrev)^2
        - CalculateDistance(xNext, x)^2
        )

    if verbose
        println("Site=$Site, x=$(round(x,digits=3)), xTest=$(round(xTest,digits=3))")
        println("ΔS=$(round(ΔS,digits=3)), exp(-ΔS)=$(round(exp(-ΔS),digits=3)).")
    end

    if ΔS <= 0
        Config.Lattice[Site] = xTest
        AccSteps += 1
        if verbose
            println("  Accepted.")
        end

    elseif rand() < exp(-ΔS)
        Config.Lattice[Site] = xTest
        AccSteps += 1
    else
        if verbose
            println("  Not accepted.")
        end
    end

    return AccSteps
end

"""
Heatbath update.
"""
function HeatBathUpdate!(Config::Configuration, i::Int64)
    # Gaussian parameter linked to simulations parameters
    Alpha = 1/(Config.Eta * Config.SimBeta)

    # Interpolate lowest action postion
    Center = (Config.Lattice[i-1] + Config.Lattice[i+1])/2

    # Extract update from gaussian pdf with null center e unitary variance,
    # then recenter the result
    Update = mod(Center + randn()/(2*sqrt(Alpha)), 1)
    Config.Lattice[i] = Update

    return Config
end

"""
Cluster "tailor" update [Bonati, D'Elia 2018], with tolerance ε.
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

"""
Plot the path given by Config.Lattice, including "pacman effect".
"""
function PlotPath(Config; location="nw")

    Q = round(Int, CalculateQ(Config))

    p = plot(
        xlims=(1,Config.N),
        ylims=(0.0,1.0),
        xlabel=L"i",
        ylabel=L"x_i",
        title=L"""
            $\tilde \beta = %$(round(Config.SimBeta, digits=2)),
            N = %$(Config.N), \eta = %$(round(Config.Eta, digits=2))$
            """
    )

    PacmanUpIndices = findall(diff(Config.Lattice) .< -0.5)
    PacmanDownIndices = findall(diff(Config.Lattice) .> 0.5)

    for i in 1:length(Config.Lattice)-1
        if i in PacmanUpIndices
            plot!([i, i+1], [Config.Lattice[i], Config.Lattice[i+1]+1], color=MyColors[4], label=nothing)
            plot!([i, i+1], [Config.Lattice[i]-1, Config.Lattice[i+1]], color=MyColors[4], label=nothing)
        elseif i in PacmanDownIndices
            plot!([i, i+1], [Config.Lattice[i], Config.Lattice[i+1]-1], color=MyColors[1], label=nothing)
            plot!([i, i+1], [Config.Lattice[i]+1, Config.Lattice[i+1]], color=MyColors[1], label=nothing)
        else
            plot!([i, i+1], [Config.Lattice[i], Config.Lattice[i+1]], color="black", label=nothing)
        end
    end

    if location=="nw"
        annotate!((0.02,0.92), text(L"$Q=%$Q$", :left, 10))
    elseif location=="ne"
        annotate!((0.98,0.92), text(L"$Q=%$Q$", :right, 10))
    end

    return p
end

"""
Plot the tailor update as a dashed line with highlighted endpoints. Write over current plot.
xxNewPlot is the output of the TailorUpdate! function.
"""
function PlotTailorUpdate!(Config, i, iEnd, xxNewPlot)
    hline!([Config.Lattice[i]], color="gray", alpha=0.3, label=nothing)

    PacmanUpIndices = findall(diff(xxNewPlot) .< -0.5)
    PacmanDownIndices = findall(diff(xxNewPlot) .> 0.5)

    for j in i:iEnd # routine identical to PlotPath
         index = j-i+1 # varies between 1:(iEnd-i+1)
         if index in PacmanUpIndices
             plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]+1], color="black", label=nothing, linestyle=:dash)
             plot!([j, j+1], [xxNewPlot[index]-1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
        elseif index in PacmanDownIndices
            plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]-1], color="black", label=nothing, linestyle=:dash)
            plot!([j, j+1], [xxNewPlot[index]+1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
        else
            plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]], color="black", linestyle=:dash, label=nothing)
        end
    end

    scatter!([i, iEnd+1], [Config.Lattice[i], Config.Lattice[iEnd+1]], color="black",
        markersize=2, label=nothing)
end


function main()
    NMetro = 10000
    N = 30

    Config = SetLattice(1.0, N)
    Counter = 0

    println("Performing $NMetro Metropolis steps...")

    for i in 1:NMetro
        Site = i % (N-2) + 2 # from 2 to N-1 (sequential)
        # Site = rand(2:N-1) # (random)
        Counter = MetropolisUpdate!(Config, Site, Counter; Δ=0.05, verbose=false)
    end

    println("Accepted steps: $Counter/$NMetro")

    p = PlotPath(Config; location="nw")

    println("Performing tailor update...")
    i = 2
    ε = 0.1
    iEnd, xxNewPlot = TailorUpdate!(Config, i, ε)

    PlotTailorUpdate!(Config, i, iEnd, xxNewPlot)

end

main()
