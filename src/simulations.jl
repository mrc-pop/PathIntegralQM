#!/usr/bin/julia

using DelimitedFiles
using Dates
using UnicodePlots

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/processing.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

"""
Simulation file. The settings are imported from /setup/simulations_setup.jl
"""
function main()

    println()

    for (SizeIndex, N) in enumerate(NN), SimBeta in SimBetas

        TailorStep = TailorSteps[SizeIndex]

        Scheme = Heatbath ? "Heatbath" : "Metropolis"

        println()

        @info "\nModel settings" N SimBeta
        @info Scheme*" settings" Heatbath Δ Sequential NSweepsTherm NSweeps
        @info "Tailor settings" TailorStep ε_over_η

        # Prepare file and write header
        DirPathOut = PROJECT_ROOT * "/../simulations/N=$N/"
        FilePathOut = DirPathOut * "SimBeta=$SimBeta.txt"
        mkpath(dirname(FilePathOut))
        open(FilePathOut, "w") do io
            write(io, "# Δ=$Δ, Sequential=$Sequential, NSweepsTherm=$NSweepsTherm\n")
            write(io, "# TailorStep=$TailorStep, ε_over_η=$ε_over_η\n")
            write(io, "# [Calculated at $(now())]\n")
        end

        # Initialize lattice and counter
        Config = SetLattice(SimBeta, N)
        ε = ε_over_η * Config.Eta
        Counter = 0
        CounterTailor = [0, 0, 0]

        # Thermalization
        println("\nPerforming $NSweepsTherm " * Scheme * " sweeps for thermalization...")
        @time for i in 1:NSweepsTherm
            CurrentSweepSteps = N*(i-1)
            for j in 1:N
                
                # Choose site
                if Sequential
                    Site = mod1(CurrentSweepSteps + j, N)
                elseif !Sequential
                    Site = rand(1:N)
                end
                
                # Perform update
                if Heatbath
                    HeatBathUpdate!(Config, Site)
                elseif !Heatbath
                   MetropolisUpdate!(Config, Site; Δ) 
                end
                
                # Perform tailor update
                if TailorStep !== 0 && mod(CurrentSweepSteps + j, TailorStep) == 0
                    FoundT, AccT, _ = TailorUpdate!(Config, rand(1:N), ε)
                    CounterTailor += [FoundT, AccT, 1]
                end
            end
        end

        unicodeplots()
        p = PlotPathUnicode(Config)
        @info "Plot after thermalization " p

        #gui()

        # Local update sweeps
        println("\nPerforming $NSweeps "* Scheme * " sweeps of the whole lattice...")

        NMetro = NSweeps * N
        QQ = zeros(NSweeps)

        @time for i in 1:NSweeps

            CurrentSweepSteps = N*(i-1)

            # Sweep over lattice
            for j in 1:N

                if Sequential
                    Site = mod1(CurrentSweepSteps + j, N)
                elseif !Sequential
                    Site = rand(1:N)
                end

                # Perform update
                if !Heatbath
                    Accepted = MetropolisUpdate!(Config, Site; Δ)
                    Counter += Accepted
                elseif Heatbath
                    HeatBathUpdate!(Config, Site)
                    Counter += 1
                end

                # Every TailorStep updates, perform a tailor update.
                # if rand(1:N) == 1 # TODO change
                if TailorStep !== 0 && mod(CurrentSweepSteps + j, TailorStep) == 0
                    Site = rand(1:N)
                    FoundT, AccT, _ = TailorUpdate!(Config, Site, ε)
                    CounterTailor += [FoundT, AccT, 1]
                end
            end

            # Save the value of Q after the sweep
            QQ[i] = CalculateQ(Config)
        end

        # Save QQ on file by appending
        open(FilePathOut, "a") do io
            writedlm(io, QQ, '\n')
        end

        println("\nSaved QQ data to $FilePathOut.")

        println("Local update acceptance: $(round(Counter/NMetro, digits=3)).")
        FoundTRatio = round(CounterTailor[1] / CounterTailor[3], digits=3)
        AccTRatio = round(CounterTailor[2] / CounterTailor[1], digits=3)
        println("Number of tailor updates: $(CounterTailor[3]).")
        println("Fraction of times iEnd was found: $FoundTRatio; "*
            "out of these, acceptance was $AccTRatio.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
