PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")

const N = 10
const SimBeta = 2.0
const Δ = 0.5

function main()

    Config = SetLattice(SimBeta, N)

    for Site in 1:N
        MetropolisUpdate!(Config, Site)
    end

    PlotPath(Config)

    Site = N-1
    Found, Acc, iEnd, xxNewPlot = TailorUpdate!(Config, Site, SimBeta/N)
    println("i = $Site, iEnd, = $iEnd")
    println("xxNewPlot = $xxNewPlot")
    println("Lattice = $(Config.Lattice)")

    plot!()

    iEnd !== nothing && PlotTailorUpdate!(Config, Site, iEnd, xxNewPlot)
end

main()
