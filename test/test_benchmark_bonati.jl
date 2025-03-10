#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

"""
Code to reproduce prof. Bonati's benchmark for N=10, N=12, and check if
our metropolis results and times are consistent.
"""

const SimBeta = 2.0
const NN = [10, 12]

const NSweepsTherm = Int(1e2)
const NSweeps = Int(1e8)
const Δ = 0.5
const Sequential = true

const TailorStep = 20
const ε = 0.01

const k = 10000

function main()

    for N in NN
        printstyled("\nWorking on N=$N", color=:yellow)
        println("\nAlgorithm settings: Δ=$Δ, Sequential=$Sequential.")
        TailorStep !== 0 && println("Tailor update every $TailorStep metro updates.")

        # Initialize lattice and counters
        Config = SetLattice(SimBeta, N)
        Counter = 0
        FoundTailor = 0

        # Thermalization
        println("\nPerforming $NSweepsTherm Metropolis sweeps for thermalization...")
        for i in 1:(NSweepsTherm*N)
            if Sequential
                Site = (i-1) % (N-2) + 2
            else
                Site = rand(2:N-1)
            end
            MetropolisUpdate!(Config, Site; Δ)
        end

        # Metropolis sweeps
        println("Lattice parameters: N=$N, SimBeta=$SimBeta, eta^2=$(round((SimBeta/N)^2, digits=5))")
        println("Performing $NSweeps Metropolis sweeps of the whole lattice...")

        NMetro = NSweeps * N
        QQ = fill(0, NSweeps)

        for i in 1:NSweeps

            CurrentMetroSteps = N*(i-1)

            # Sweep over lattice
            for j in 1:N

                if Sequential
                    Site = (j-1) % (N-2) + 2
                else
                    Site = rand(2:N-1)
                end

                # Perform the Metropolis update
                Accepted = MetropolisUpdate!(Config, Site; Δ)
                Counter += Accepted

                # Perform the tailor update.
                # CurrentMetroSteps + j = how many Metro steps have been already performed
                if TailorStep !== 0 && mod(CurrentMetroSteps + j, TailorStep) == 0
                    Found, _ = TailorUpdate!(Config, rand(1:N), ε)
                    FoundTailor += Found
                end
            end

            # Save the value of Q after the sweep
            QQ[i] = round(Int,CalculateQ(Config))
        end

        println("Accepted steps: $Counter/$NMetro (acceptance $(round(Counter/NMetro, digits=2))).")

        if TailorStep !== 0
            NTailor = round(Int64, NMetro/TailorStep)
            println("How many times tailor found iEnd: $FoundTailor/$NTailor "*
                "(ratio $(round(FoundTailor/NTailor, digits=2)))")
        end

        println("\nBlocking Q and Q² with k=$k, number of blocks: $(round(Int64, length(QQ)/k))")

        Q, SigmaQ, Q2, SigmaQ2 = BlockQ(QQ, k)

        println("⟨Q⟩ = $Q ± $SigmaQ")
        println("⟨Q²⟩ = $Q2 ± $SigmaQ2")
    end
end

@time main()
