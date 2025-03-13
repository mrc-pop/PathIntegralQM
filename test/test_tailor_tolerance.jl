#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")
using UnicodePlots
using DelimitedFiles

const NSweepsTherm = Int(1e3)
const NSweeps = Int(1e3)                    # number of updates of the whole lattice
const Δ = 0.5
const sequential = false # (Metropolis)
const NN = [400]                            # number of lattice points
const SimBeta = 2.0
const εε_over_ηη = [1.0]     # tolerance of tailor method, in units of η
const MeasureQEvery = 1
const TailorEvery = 10
const FilePathOut = PROJECT_ROOT * "/tailor_tolerance/NN=$(NN)_SimBeta=$SimBeta.txt"

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

            Config = SetLattice(SimBeta, N)
            FoundTailor = 0
            AccTailor = 0
            MeasureCount = 0
            TailorCount = 0

            for j in 1:(NSweepsTherm*N)
                Site = sequential ? mod1(j, N) : rand(1:N) # compact if-else
                MetropolisUpdate!(Config, Site; Δ)
            end

            LL = fill(0, NSweeps*N) # lengths
            QQ = fill(0, Int(ceil(N*NSweeps/MeasureQEvery)))

            if FilePathOut!=""
                open(FilePathOut, "a") do io
                    println("Opening file...")
                    write(io, "# NSweeps = $NSweeps")
                    write(io, "# N, SimBeta, NSweeps, ε_over_η, TailorEvery $(now())")
                end
            end

            println("Performing $(round(NSweeps,sigdigits=3)) sweeps of the lattice... ")

            @time for i in 1:NSweeps

                # Sweep over lattice using Metropolis
                for j in 1:N

                    if mod((N-1)*i+j, MeasureQEvery) == 0
                        MeasureCount += 1
                        QQ[MeasureCount] = CalculateQ(Config)
                    end

                    Site = sequential ? mod1(j, N) : rand(1:N)
                    MetropolisUpdate!(Config, Site; Δ)

                    if mod((N-1)*i+j, TailorEvery) == 0
                        TailorCount += 1
                        Site = rand(1:N)
                        Found, Acc, iEnd, _ = TailorUpdate!(Config, Site, ε; verbose=false)

                        # If iEnd is found, store length
                        if iEnd !== nothing
                            LL[TailorCount] = mod1(iEnd - Site, N)
                        end

                        # Store outcome of tailor update
                        FoundTailor += Found
                        AccTailor += Acc
                    end
                end
            end

            # exclude zeros (i.e. when iEnd was not found)
            LL = LL[LL .!= 0]

            MeanL = mean(LL)
            StdL = std(LL)

            FoundRatio = round(FoundTailor/TailorCount, digits=3)
            AccRatio = round(AccTailor/FoundTailor, digits=3)

            println("""
             Over $(round(TailorCount, sigdigits=3)) tailor updates, I found iEnd $(round(FoundTailor, sigdigits=3)) times (ratio $FoundRatio),
             and accepted the tailor metropolis step $AccTailor times (acceptance $AccRatio).
             Length of the cluster: L = $(round(MeanL, digits=2)) ± $(round(StdL,digits=2))
            """)

            # Calculate Q
            kk = 1:500:10000
            EQ2_Blocks = []
            for k in kk
                QQ = QQ[1:MeasureCount]
                _, _, _, EQ2 = BlockQ(QQ, k)
                push!(EQ2_Blocks, EQ2)
            end

            p = lineplot(kk, EQ2_Blocks, width=150, height=30)
            @info p

            print("Choose your k among $kk: ")
            Waiting = true
            Input = readline()
            SelectedK = parse(Int, Input)
            while Waiting
                if SelectedK ∈ collect(kk)
                    Waiting = false
                else
                    print("Selected k is not in the range of kk. Try again: ")
                    Input = readline()
                    SelectedK = parse(Int, Input)
                end
            end

            SelectedIndex = findfirst(SelectedK .== kk)

            Q, EQ, Q2, EQ2 = BlockQ(QQ, SelectedK)
            println("""
            For ε/η = $Fraction, k = $SelectedK, TailorEvery = $TailorEvery:
                FoundRatio = $(round(FoundRatio, sigdigits=3)),
                AccRatio = $(round(AccRatio, sigdigits=3)),
                Average Q = $Q ± $EQ,
                Average Q² = $Q2 ± $EQ2
            """)

            FilePathOut!="" && open(FilePathOut, "a") do io
                writedlm(io, "$N; $Fraction")
            end
        end
    end
end

main()
