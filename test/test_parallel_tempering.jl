#!/usr/bin/julia

"""
Code to determine the optimal value of SimBeta_max = (Ratio)^(NR-1) * SimBeta for simulating
the system with parallel tempering. Simulation settings are NOT imported from other scripts.
"""

using Printf
include(joinpath(@__DIR__, "../src/simulations_parallel_tempering.jl"))

PROJECT_ROOT = @__DIR__

const NN = [300]
const SimBetas = [2.0]
const SimBeta = SimBetas[1]

const Heatbath = false
# const NSweepsTherm = Int(1e2)
const NSweeps = Int(1e6)
const Δ = 0.5
const Sequential = true

const TailorSteps = round.(0.0 .* NN)
const ε_over_η = 0.2

const MeasureEvery = 1 .* NN

"""PT settings."""
const NR = 5
const Ratios = [1.02, 1.03, 1.04, 1.05, 1.06, 1.08, 1.10, 1.14, 1.16, 1.18] #[1.012, 1.025, 1.05, 1.075, 1.1] # CHANGE! # 1.01:0.01:1.08
const SwapStep = 600

"""Analysis settings."""
const LengthsBlockSizes = Dict(
    # N => k
    200 => 500,
    300 => 2000,
    400 => 2000
)


function main()
    if length(ARGS) != 1
        @error """
        Wrong number of arguments. How to use this script?
        julia test_parallel_tempering.jl <Mode>
        where the mode can be --simulations, --blocking, --analysis
        """
        return
    end

    Mode = ARGS[1]

    if Mode == "--simulations"
        runSimulations()
    elseif Mode == "--blocking"
        runBlocking()
    elseif Mode == "--analysis"
        runAnalysis(LengthsBlockSizes[NN[1]])
    else
        @error "Invalid argument. Use julia test_parallel_tempering.jl <Mode>."
    end
end

function runSimulations()
    for K in Ratios
        FolderOut = PROJECT_ROOT * "/test_parallel_tempering_data/NR=$(NR)_Ratio=$K/"
        mkpath(FolderOut)
        PTRoutine(NR, K, SwapStep, FolderOut)
    end
end

function runBlocking(;SavePlot=true)

    BlockSizes = vcat(1, 5:5:20, 40:20:100, 150:50:300, 400:200:4000)

    plot()
    for K in Ratios
        for N in NN
            FolderOut = PROJECT_ROOT * "/test_parallel_tempering_data/NR=$(NR)_Ratio=$K/"
            FilePathIn = FolderOut * "/N=$N/SimBeta=$SimBeta.txt"

            QQ = readdlm(FilePathIn, '\n', Int64; comments=true)
            QQ2 = QQ .^ 2
            MeanQ2 = mean(QQ2)

            @info "Blocking for different block sizes" NR K N SimBeta MeanQ2

            ErrorsQ = fill(0.0, length(BlockSizes))
            ErrorsQ2 = fill(0.0, length(BlockSizes))

            @time for (i, k) in enumerate(BlockSizes)
                @info "Blocking Q and Q² with k=$k, number of blocks: $(round(Int64, length(QQ) / k))"
                _, ErrorsQ[i], _, ErrorsQ2[i] = BlockQ(QQ, k)
            end

            plot!(BlockSizes, ErrorsQ2, xlabel=L"k", ylabel=L"\sigma_{Q^2}",
                title=L"$\tilde\beta = %$SimBeta, N_R = %$NR$ ", label=L"K=%$K", legend=:bottomright)
        end
    end

    if SavePlot
        savefig(PROJECT_ROOT * "/test_parallel_tempering_data/blocking/blocking_NN=$(NN)_SimBeta=$(SimBeta)_NR=$(NR).pdf")
        println("Saved figure to /test_parallel_tempering_data/blocking/blocking_NN=$(NN)_SimBeta=$(SimBeta)_NR=$(NR).pdf")
    end
end

function runAnalysis(k)

    FilePathOut = PROJECT_ROOT * "/test_parallel_tempering_data/SimBeta=$(SimBeta)_analysis.txt"
    open(FilePathOut, "a") do io
        write(io, "# NR=$NR, Block size=$k, NSweeps=$NSweeps, SwapStep=$SwapStep [calculated $(now())]\n")
        write(io, "# N, K, ErrorQ2\n")
    end

    ErrorsQ2 = Float64[]

    for K in Ratios
        for N in NN
            @info "Blocking N=$N with block size $k, Ratio=$K..."
            Folder = PROJECT_ROOT * "/test_parallel_tempering_data/NR=$(NR)_Ratio=$K/"
            FilePathIn = Folder * "/N=$N/SimBeta=$SimBeta.txt"

            QQ = readdlm(FilePathIn, '\n', Int64; comments=true)

            Q, ErrorQ, Q2, ErrorQ2 = BlockQ(QQ, k)

            # Write on file NN, ErrorQ2, K
            open(FilePathOut, "a") do io
                write(io, "$N, $K, $ErrorQ2\n")
            end

            push!(ErrorsQ2, ErrorQ2*sqrt(NR))
        end
    end

    plot(Ratios, ErrorsQ2, xlabel=L"K", ylabel=L"\sqrt{N_R} \cdot \overline\sigma_{\overline\chi}", label=L"$N_R = %$NR$",
        title=L"Search for optimal $K$", legend=:bottomright)

    savefig(PROJECT_ROOT * "/test_parallel_tempering_data/result_NN=$(NN)_SimBeta=$(SimBeta)_NR=$(NR).pdf")
    println("Saved figure to /test_parallel_tempering_data/result_NN=$(NN)_SimBeta=$(SimBeta)_NR=$(NR).pdf")
end

@time main()
