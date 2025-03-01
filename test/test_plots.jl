#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")

const N = 50
const SimBeta = 2.0
const Δ = 0.2

function main()

    println("""
    Test of Metropolis updates and plot of the resulting path; tailor update and plot
    of the resulting path, with an indication of whether it has been accepted.
    """)

    Config = SetLattice(SimBeta, N)

    for i in 1:(10*N)
        Site = mod1(i,N)
        MetropolisUpdate!(Config, Site; Δ)
    end

    println("Before tailor, Q=$(CalculateQ(Config))")

    p = PlotPath(Config)

    savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$N.pdf")

    plot(p)

    Site = 10
    Found, Acc, iEnd, xxNewPlot = TailorUpdate!(Config, Site, 0.02)

    println("iEnd=$iEnd. Found? $(Bool(Found)). Accepted? $(Bool(Acc))")

    PlotTailorUpdate!(Config, Site, iEnd, xxNewPlot)

    println("After tailor, Q=$(CalculateQ(Config))")

    if Bool(Found) == 1
        savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$(N)_tailor.pdf")
        PlotPath(Config)
        savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$(N)_after.pdf")
    end

end

main()
