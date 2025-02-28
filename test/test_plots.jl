#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")

const N = 50
const SimBeta = 2.0
const Δ = 0.5

function main()

    println("""
    Test of Metropolis updates and plot of the resulting path; tailor update and plot
    of the resulting path, with an indication of whether it has been accepted.
    """)

    Config = SetLattice(SimBeta, N)

    for i in 1:(100*N)
        Site = mod1(i,N)
        MetropolisUpdate!(Config, Site; Δ)
    end

    println("Before tailor, Q=$(CalculateQ(Config))")

    PlotPath(Config)

    Site = N-1
    Found, Acc, iEnd, xxNewPlot = TailorUpdate!(Config, Site, 0.2*SimBeta/N)

    println("iEnd=$iEnd. Found? $(Bool(Found)). Accepted? $(Bool(Acc))")

    PlotTailorUpdate!(Config, Site, iEnd, xxNewPlot)

    println("After tailor, Q=$(CalculateQ(Config))")

    gui()

    plot!()
end

main()
