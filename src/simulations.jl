#!/usr/bin/julia

using DelimitedFiles
using Dates
using UnicodePlots

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

    for (sizeindex, N) in enumerate(NN), SimBeta in SimBetas

        TailorStep = TailorSteps[sizeindex]

        AlgoName = Heatbath ? "Heatbath" : "Metropolis"

        println()

        @info "\nModel settings" N SimBeta
#        @info AlgoName*" settings" Heatbath Δ Sequential NSweepsTherm NSweeps
        @info "Tailor settings" TailorStep ε_over_η

        # Prepare file and write header
        FolderOut = PROJECT_ROOT * "/../simulations/N=$N/"
        FilePathOut = FolderOut * "SimBeta=$SimBeta.txt"
        mkpath(dirname(FilePathOut))
        open(FilePathOut, "w") do io
            write(io, "# Δ=$Δ Sequential=$Sequential NSweepsTherm=$NSweepsTherm\n")
            write(io, "# TailorStep=$TailorStep ε_over_η=$ε_over_η\n")
            write(io, "# [Calculated at $(now())]\n")
        end

        # Initialize lattice and counter
        Config = SetLattice(SimBeta, N)
        ε = ε_over_η * Config.Eta
        Counter = 0
        CounterTailor = [0, 0, 0]

        @info AlgoName*" settings" Heatbath Δ Sequential NSweepsTherm NSweeps

        # Thermalization
        println("\nPerforming $NSweepsTherm " * AlgoName * " sweeps for thermalization...")
        @time for i in 1:NSweepsTherm
            CurrentMetroSteps = N*(i-1)
            for j in 1:N
                # Choose site
                if Sequential
                    Site = mod1(CurrentMetroSteps + j, N)
                else
                    Site = rand(1:N)
                end
                # Perform update
                if Heatbath == false
                    MetropolisUpdate!(Config, Site; Δ)
                else
                    HeatBathUpdate!(Config, Site)
                end
                # Perform tailor update
                if TailorStep !== 0 && mod(CurrentMetroSteps + j, TailorStep) == 0
                    TailorUpdate!(Config, rand(1:N), ε)
                end
            end
        end

        unicodeplots()
        p = PlotPathUnicode(Config)
        @info "Plot after thermalization " p

        # Local update sweeps
        println("\nPerforming $NSweeps "* AlgoName * " sweeps of the whole lattice...")

        NMetro = NSweeps * N
        MeasureInterval = QSteps[sizeindex]
        ExpectedMeasurements = NSweeps * N ÷ MeasureInterval
        QQ = fill(0, ExpectedMeasurements)
        MeasureCount = 0

        Site = 1 # starting site

        @time for i in 1:NSweeps

            CurrentMetroSteps = N*(i-1)

            # Sweep over lattice
            for j in 1:N

                if Sequential == false
                    Site = rand(1:N)
                end

                # Perform update
                if Heatbath == false
                    Accepted = MetropolisUpdate!(Config, Site)
                    Counter += Accepted
                    Site = mod1(Site+1, N)
                else
                    HeatBathUpdate!(Config, Site)
                    Counter += 1
                end

                #@info "Performed local update." PlotPathUnicode(Config)

                # Every TailorStep updates, perform a tailor update.
                if TailorStep !== 0 && mod(CurrentMetroSteps + j, TailorStep) == 0
                    FoundT, AccT, iEnd, _ = TailorUpdate!(Config, Site, ε)
                    CounterTailor += [FoundT, AccT, 1]
                    #@info "Performed tailor update. Found? $FoundT. i=$Site, iEnd=$iEnd" PlotPathUnicode(Config)

                    if iEnd !== nothing
                        Site = mod1(iEnd + 1, N)
                    end
                end

                # Measure Q
                if mod(CurrentMetroSteps + j, MeasureInterval) == 0
                    MeasureCount += 1
                    QQ[MeasureCount] = CalculateQ(Config)
                end

            end
        end

        # Ensure we only save the measured values
        QQ = QQ[1:MeasureCount]

        printstyled("\nAverage of Q² over $MeasureCount measurements: $(mean(QQ.^2))\n", color=:yellow)

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
