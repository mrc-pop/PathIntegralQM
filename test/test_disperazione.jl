PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

using UnicodePlots

"""
Trying to reproduce `Q vs time` plot of Prof. Bonati.
"""

# Model settings
const SimBeta = 2.0
const N = 20
# const QForce = 0 # Force the initial lattice to have winding number QForce

# MCMC update settings
const NSweepsTherm = Int(1e4)
const NSweeps = 10
const heatbath = true
const Δ = 0.5
const TailorEvery = 0 # [0 = none]

# Measurement settings
const MeasureEvery = 1

function main()

    Config = SetLattice(SimBeta, N)

    # if QForce == 0
    #     # Set random Q=0 path
    #     Config.Lattice .= [0.5*rand() + 0.25 for i in 1:N]
    # else
    #     # Set uniform motion path
    #     Config.Lattice .= [mod( 0.5 + j * QForce/N ,1) for j in 1:N]
    # end

    Config.Lattice .= [rand() for _ in 1:N]

    @info "Model parameters" SimBeta N Config.Eta
    @info "Initial config" PlotPathUnicode(Config)

    Counter = 0
    QQ = Int64[]

    for i in 1:NSweepsTherm
        for j in 1:N
            MetropolisUpdate!(Config, rand(1:N); Δ)
        end
    end

    p = PlotPathUnicode(Config)

    @info "Config after thermalization" PlotPathUnicode(Config)

    for i in 1:NSweeps

        CurrentMetroSteps = N*(i-1)

        # Sweep over lattice
        for j in 1:N

            Site = rand(1:N)

            # Perform update
            if heatbath == false
                Accepted = MetropolisUpdate!(Config, Site; Δ)
                Counter += Accepted
            else
                HeatBathUpdate!(Config, Site)
                Counter += 1
            end

            # Tailor update
            if TailorEvery !== 0 && mod(CurrentMetroSteps + j, TailorEvery) == 0
                Site = rand(1:N)
                Found, Acc, _ = TailorUpdate!(Config, Site, 0.2 * Config.Eta)
            #    println("On $i-th sweep, j=$j: Q=$(CalculateQ(Config)) -> tailor: Found=$Found, Acc=$Acc")
            end

            if mod(CurrentMetroSteps + j, MeasureEvery) == 0
                Q = CalculateQ(Config)
                push!(QQ, Q)
            end
        end
    end

    q = lineplot(MeasureEvery.* (1:length(QQ)), QQ, xlabel="time", ylabel="Q", width=:120, height=:30)

    @info "Q vs Metropolis step (one measurement every $MeasureEvery local updates)" q mean(QQ) var(QQ) mean(QQ.^2)
end

@time main()
