#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

const NSweepsTherm = Int(1e3)
const NSweeps = Int(1e4)                    # number of updates of the whole lattice
const Δ = 0.3
const sequential = true
const NN = [100, 200, 300, 400]                            # number of lattice points
const SimBeta = 2.0
const εε_over_ηη = [0.2, 1.0, 2.0, 10.0]     # tolerance of tailor method, in units of η

function main()

    println("""
    \nPerforming test of tailor update success of finding iEnd and acceptance,
    for different values of the tolerance ε, in units of η. (sequential = $sequential)
    """)

    for N in NN
        printstyled("\nWorking on N=$N", color=:yellow)

        η = SimBeta/N

        for Fraction in εε_over_ηη
            ε = Fraction * η
            println("\nUsing N=$N, ε/η = $Fraction (ε=$ε, η=$η)")

            Config = SetLattice(1.0, N)
            FoundTailor = 0
            AccTailor = 0

            for j in 1:(NSweepsTherm*N)
                Site = sequential ? mod1(j, N) : rand(1:N) # compact if-else
                MetropolisUpdate!(Config, Site; Δ)
            end

            LL = fill(0, NSweeps)

            for i in 1:NSweeps

                # Sweep over lattice using Metropolis
                for j in 1:N
                    Site = sequential ? mod1(j, N) : rand(1:N)
                    MetropolisUpdate!(Config, Site; Δ)
                end

                # Perform tailor on a random site
                Site = rand(1:N)
                Found, Acc, iEnd, _ = TailorUpdate!(Config, Site, ε; verbose=false)

                # If iEnd is found, store length
                if iEnd !== nothing
                    LL[i] = mod1(iEnd - Site, N)
                end

                # Store outcome of tailor update
                FoundTailor += Found
                AccTailor += Acc
            end

            # exclude zeros (i.e. when iEnd was not found)
            LL = LL[LL .!= 0]

            MeanL = mean(LL)
            StdL = std(LL)

            FoundRatio = round(FoundTailor/NSweeps, digits=3)
            AccRatio = round(AccTailor/FoundTailor, digits=3)

            println("""
            Over $NSweeps tailor updates, I found iEnd $FoundTailor times (ratio $FoundRatio),
            and accepted the tailor metropolis step $AccTailor times (acceptance $AccRatio).
            Length of the cluster: L = $(round(MeanL, digits=2)) ± $(round(StdL,digits=2))
            """)
        end
    end
end

@time main()
