#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")
using UnicodePlots
using DelimitedFiles
using Dates

const NSweepsTherm = Int(1e3)
const NMetro = Int(400 * 1e5)                    # number of updates of the whole lattice
const Δ = 0.5
const sequential = false # (Metropolis)
const NN = [400]#[75, 100, 125, 150, 175, 200, 300, 400]  # number of lattice points
const SimBeta = 2.0
const εε_over_ηη = [1.0]     # tolerance of tailor method, in units of η
const MeasureQEvery = 1
const TailorEverys = [50]
const FilePathOut = PROJECT_ROOT * "/tailor_tolerance/simulations/NMetro=$(round(NMetro,sigdigits=2))_SimBeta=$(SimBeta).txt"

function TailorSimulation()

    println("""
    \nPerforming test of tailor update success of finding iEnd and acceptance,
    for different values of the tolerance ε, in units of η. (sequential = $sequential)
    """)

    if FilePathOut!=""
        open(FilePathOut, "a") do io
            println("Opening file...")
            write(io, "# SimBeta=$SimBeta, sequential=$sequential [calculated $(now())]\n")
            write(io, "# N; NSweeps; ε_over_η; k (BlockSize); TailorEvery; MeasureQEvery; ⟨Q⟩; e⟨Q⟩; ⟨Q²⟩; e⟨Q²⟩; Δt_Metro [s]; Δt_Blocking [s]\n")
        end
    end

    for N in NN, TailorEvery in TailorEverys
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

            LL = fill(0, NMetro) # lengths
            QQ = fill(0, Int(ceil(NMetro/MeasureQEvery)))

            NSweeps = floor(Int, NMetro/N)
            println("Performing $NSweeps sweeps of the lattice... ")

            ElapsedTimeSweeps = @elapsed begin
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
            println("Calculating std dev for different block sizes...")
            ElapsedTimeBlocking = @elapsed begin
            kk = 1:500:10000
            EQ2_Blocks = []
            @time for k in kk
                QQ = QQ[1:MeasureCount]
                _, _, _, EQ2 = BlockQ(QQ, k)
                push!(EQ2_Blocks, EQ2)
            end
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
                write(io, "$N; $NSweeps; $Fraction; $SelectedK; $TailorEvery; $MeasureQEvery; $Q; $EQ; $Q2;",
                " $EQ2; $ElapsedTimeSweeps; $ElapsedTimeBlocking\n")
            end
        end
    end
end

function TailorSimulationPlot(FilePath::String)
    Data = readdlm(FilePath, ';', comments=true)

    Sizes = Data[:,1]
    Q = Data[:,6]
    EQ = Data[:,7]
    Q2 = Data[:,8]
    EQ2 = Data[:,9]

    @info "Plot parameters" NN SimBeta

    # Create plot
    ExponentTitle = round(Int,log10(NMetro))
    p = scatter(Sizes, Q2,
            title=L"\tilde\beta = %$SimBeta, N_\mathrm{metro} = 10^{%$ExponentTitle}, n_\mathrm{tailor}=%$TailorEvery",
            markersize=2,
            yerror=EQ2,
            xlabel=L"N",
            ylabel=L"\langle Q^2 \rangle",
            label=:none
        )

    pgfplotsx()

    # Add theoretical line
    if SimBeta > 1.0
        Plots.hline!([SimBeta],
            label=:none,#L"$\tilde\beta$",
            color="gray",
            linestyle=:dash)
    elseif SimBeta < 0.5
        Plots.hline!([exp(-1 / (2*SimBeta))],
            label=:none,#L"$\exp(-1/2\tilde\beta)$",
            color="gray",
            linestyle=:dash)
    end

    # Save plot
    Plots.savefig(p, PROJECT_ROOT * "/tailor_tolerance/plots/Q2_vs_N_NMetro=$(NMetro)_SimBeta=$(SimBeta).pdf")
    println("Plot of Q2 saved")

    # Create plot of Q vs NN
    p2 = scatter(Sizes, Q,
            title=L"\tilde\beta = %$SimBeta, N_\mathrm{metro} = 10^{%$ExponentTitle}, n_\mathrm{tailor}=%$TailorEvery",
            markersize=2,
            yerror=EQ,
            xlabel=L"N",
            ylabel=L"\langle Q \rangle",
            label=:none)
    Plots.hline!([0.0],
            label=:none,
            color="gray",
            linestyle=:dash)

    # Save plot
    Plots.savefig(p2, PROJECT_ROOT * "/tailor_tolerance/plots/Q_vs_N_NMetro=$(NMetro)_SimBeta=$(SimBeta).pdf")
    println("Plot of Q saved")
end

function main()
    print("Simulation (s) or plot (p)?: ")
    Waiting = true
    Input = readline()
    while Waiting
        if Input=="s"
            TailorSimulation()
            Waiting = false
        elseif Input=="p"
            TailorSimulationPlot(FilePathOut)
            Waiting = false
        end
    end
end

main()
