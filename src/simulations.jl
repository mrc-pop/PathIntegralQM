#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/processing.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

const NSweepsTherm = 100
const NSweeps = Int(1e8)       # number of updates of the whole lattice
const Δ = 0.3
const sequential = false
const N = 100                    # number of lattice points
const SimBeta = 2.0
const ε = 0.2 * (SimBeta/N)      # tolerance of tailor method (article)
const k = Int(1e4)               # block length

function main()

    println("\nAlgorithm settings: Δ=$Δ, sequential=$sequential, ε=$ε.")

    # Initialize lattice and counter
    Config = SetLattice(SimBeta, N)
    Counter = 0
    CounterTailor = [0, 0]

    # Thermalization
    println("\nPerforming $NSweepsTherm Metropolis sweeps for thermalization...")
    for i in 1:(NSweepsTherm*N)
        if sequential
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

        # Sweep over lattice
        for j in 1:N

            if sequential
                Site = (j-1) % (N-2) + 2
            else
                Site = rand(2:N-1)
            end

            # Perform the Metropolis update
            Accepted = MetropolisUpdate!(Config, Site; Δ)
            Counter += Accepted

            # Every XXX updates, perform a tailor update
            if j % 50 == 0
                AccT, _ = TailorUpdate!(Config, Site, ε)
                CounterTailor += [AccT, 1]
            end
        end

        # Save the value of Q after the sweep
        QQ[i] = round(Int,CalculateQ(Config))
    end

    println("Accepted steps: $Counter/$NMetro (acceptance $(round(Counter/NMetro, digits=2))).")
    println("Number of tailor updates: $CounterTailor.")

    Q = mean(QQ)
    Q2 = mean(QQ.^2)

    println("\nBlocking Q and Q² with k=$k, number of blocks: $(round(Int64, length(QQ)/k))")

    BlockedQQ = BlockData(Float64.(QQ), k)
    SigmaQ = std(BlockedQQ, corrected=true) / sqrt(length(BlockedQQ))

    BlockedQQ2 = BlockData(Float64.(QQ.^2), k)
    SigmaQ2 = std(BlockedQQ2, corrected=true) / sqrt(length(BlockedQQ2))

    println("⟨Q⟩ = $Q ± $SigmaQ")
    println("⟨Q²⟩ = $Q2 ± $SigmaQ2")

    PlotPath(Config)
    plot!(size=(800, 400))
    gui()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
