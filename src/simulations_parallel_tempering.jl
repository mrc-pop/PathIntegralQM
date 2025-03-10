#!/usr/bin/julia

using DelimitedFiles
using Dates
# using UnicodePlots

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/processing.jl")
include(PROJECT_ROOT * "/modules/plots.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

const NR = 1 # number of systems to be simulated in parallel
const Ratio = 1.4 # ratio between the η of successive systems.
const SwapStep = 20 # how often to propose an exchange

"""
Simulation file with parallel tempering. The settings are imported from /setup/simulations_setup.jl
"""

function PTRoutine(NR, Ratio, SwapStep, FolderOut)

    println()

    for (sizeindex, N) in enumerate(NN), SimBeta in SimBetas

        TailorStep = TailorSteps[sizeindex]

        println()

        @info "\nModel settings" N SimBeta
        AlgoName = Heatbath ? "Heatbath" : "Metropolis"
        @info AlgoName*" settings" Heatbath Δ Sequential NSweepsTherm NSweeps
        @info "Tailor settings" TailorStep ε_over_η

        # Prepare file and write header
        FilePathOut = FolderOut * "N=$N/SimBeta=$SimBeta.txt"
        mkpath(dirname(FilePathOut))
        open(FilePathOut, "w") do io
            write(io, "# Δ=$Δ Sequential=$Sequential NSweepsTherm=$NSweepsTherm\n")
            write(io, "# TailorStep=$TailorStep ε_over_η=$ε_over_η\n")
            write(io, "# Parallel tempering with N_R = $NR\n")
            write(io, "# [Calculated at $(now())]\n")
        end

        # Initialize NR lattices and counter
        SimBetasParallel = [Ratio^k * SimBeta for k in 0:NR-1]
        EtasParallel = SimBetasParallel ./ N
        Configs = [SetLattice(β, N) for β in SimBetasParallel]
        εParallel = ε_over_η .* EtasParallel
        CountersParallel = [0 for _ in 1:NR]
        CountersTailor = [[0, 0, 0] for _ in 1:NR]

        @info "Parallel tempering settings" NR Ratio SimBetasParallel[end] SwapStep

        # TODO bring back thermalization
        # println("\n Instead of thermalization, I am using random initial lattices!")
        # for r in 1:NR
        #     Config = Configs[r]
        #     Config.Lattice .= [rand() for _ in 1:N]
        # end
        # println("Initial Qs, $(CalculateQ.(Configs))")

        # @info PlotPathUnicode(Configs[1])

        # Local update sweeps
        println("\nPerforming $NSweeps "* AlgoName * " sweeps of the whole lattice...")

        NMetro = NSweeps * N
        MeasureInterval = MeasureEvery[sizeindex]
        ExpectedMeasurements = NSweeps * N ÷ MeasureInterval
        QQ = fill(0, ExpectedMeasurements)
        MeasureCount = 0

        SitesParallel = [1 for _ in 1:NR] # starting site

        @time for i in 1:NSweeps

            CurrentMetroSteps = N*(i-1)

            # Sweep over lattice
            for j in 1:N

                # Single update of all systems in parallel
                for (r, Config) in enumerate(Configs)

                    Site = SitesParallel[r]

                    if Sequential == false
                        Site = rand(1:N)
                    end

                    # Perform update
                    if Heatbath == false
                        Accepted = MetropolisUpdate!(Config, Site; Δ=Δ)
                        CountersParallel[r] += Accepted
                        SitesParallel[r] = mod1(SitesParallel[r]+1, N)
                    else
                        HeatBathUpdate!(Config, Site)
                        CountersParallel[r] += 1
                    end

                    # Every TailorStep updates, perform a tailor update.
                    if TailorStep !== 0 && mod(CurrentMetroSteps + j, TailorStep) == 0
                        FoundT, AccT, iEnd, _ = TailorUpdate!(Config, Site, εParallel[r])
                        CountersTailor[r] += [FoundT, AccT, 1]
                        #@info "Performed tailor update. Found? $FoundT. i=$Site, iEnd=$iEnd" PlotPathUnicode(Config)

                        if iEnd !== nothing
                            SitesParallel[r] = mod1(iEnd + 1, N)
                        end
                    end
                end

                # Every SwapStep updates, perform an exchange.
                if mod(CurrentMetroSteps + j, SwapStep) == 0
                    if rand() < 0.5
                        for r in 1:(NR-1)
                            #println("\nProposing update between $r and $(r+1)...")
                            Acc = ParallelTemperingUpdate!(Configs[r], Configs[r+1]; verbose=false)
                        end

                    else
                        for r in NR:-1:2 # go down from NR to 2
                            #println("\nProposing update between $r and $(r-1)...")
                            Acc = ParallelTemperingUpdate!(Configs[r], Configs[r-1]; verbose=false)
                        end
                    end
                end

                # Measure Q
                if mod(CurrentMetroSteps + j, MeasureInterval) == 0
                    MeasureCount += 1
                    QQ[MeasureCount] = CalculateQ(Configs[1])
                    # since the first one is the physical one
                end

            end
        end

        # Ensure we only save the measured values
        QQ = QQ[1:MeasureCount]

        printstyled("\nAverage of Q²: $(mean(QQ.^2))\n", color=:yellow)

        # Save QQ on file by appending
        open(FilePathOut, "a") do io
            writedlm(io, QQ, '\n')
        end

        println("\nSaved QQ data to $FilePathOut.")

        for r in 1:NR
            println("\nSome statistics for r=$r.")
            Counter = CountersParallel[r]
            CounterTailor = CountersTailor[r]
            println("  Local update acceptance: $(round(Counter/NMetro, digits=3)).")
            FoundTRatio = round(CounterTailor[1] / CounterTailor[3], digits=3)
            AccTRatio = round(CounterTailor[2] / CounterTailor[1], digits=3)
            println("  Number of tailor updates: $(CounterTailor[3]).")
            println("  Fraction of times iEnd was found: $FoundTRatio; "*
                "out of these, acceptance was $AccTRatio.")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    PTRoutine(NR, Ratio, SwapStep, PROJECT_ROOT * "/../simulations/")
end
