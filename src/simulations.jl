#!/usr/bin/julia

using DelimitedFiles
using Dates

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/processing.jl")
include(PROJECT_ROOT * "/modules/plots.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

"""
Simulation file. The settings are imported from /setup/simulations_setup.jl
"""

function main()

    println()

    for N in NN, SimBeta in SimBetas

        @info "\nModel settings" N SimBeta
        @info "Metropolis settings" Δ sequential NSweepsTherm NSweeps
        @info "Tailor settings" TailorStep ε_over_η

        # Prepare file and write header
        FolderOut = PROJECT_ROOT * "/../simulations/N=$N/"
        FilePathOut = FolderOut * "SimBeta=$SimBeta.txt"
        mkpath(dirname(FilePathOut))
        open(FilePathOut, "w") do io
            write(io, "# Δ=$Δ sequential=$sequential NSweepsTherm=$NSweepsTherm\n")
            write(io, "# TailorStep=$TailorStep ε_over_η=$ε_over_η\n")
            write(io, "# [Calculated at $(now())]\n")
        end

        # Initialize lattice and counter
        Config = SetLattice(SimBeta, N)
        ε = ε_over_η * Config.Eta
        Counter = 0
        CounterTailor = [0, 0]

        # Thermalization
        println("\nPerforming $NSweepsTherm Metropolis sweeps for thermalization...")
        @time for i in 1:(NSweepsTherm*N)
            if sequential
                Site = (i-1) % (N-2) + 2
            else
                Site = rand(2:N-1)
            end
            MetropolisUpdate!(Config, Site; Δ)
        end

        # Metropolis sweeps
        println("\nPerforming $NSweeps Metropolis sweeps of the whole lattice...")

        NMetro = NSweeps * N
        QQ = fill(0, NSweeps)

        @time for i in 1:NSweeps

            CurrentMetroSteps = N*(i-1)

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

                # Every TailorStep updates, perform a tailor update.
                if TailorStep !== 0 && mod(CurrentMetroSteps + j, TailorStep) == 0
                    AccT, _ = TailorUpdate!(Config, rand(1:N), ε)
                    CounterTailor += [AccT, 1]
                end
            end

            # Save the value of Q after the sweep
            QQ[i] = round(Int,CalculateQ(Config))
        end

        # Save QQ on file by appending
        open(FilePathOut, "a") do io
            writedlm(io, QQ, '\n')
        end

        println("\nSaved QQ data to $FilePathOut.")

        println("Accepted steps: $Counter/$NMetro (acceptance $(round(Counter/NMetro, digits=2))).")
        println("Number of tailor updates: $CounterTailor.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
