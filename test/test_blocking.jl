#!/usr/bin/julia

using DelimitedFiles

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

"""
Code to find the optimal block length.
"""

# Import variables from setup file
include(PROJECT_ROOT * "/../src/setup/simulations_setup.jl")

# Define custom block sizes variable
const BlockSizes = vcat(1, 5:5:20, 40:20:100, 150:50:300, 400:100:1000)
const NNBlock = NN
# [1, 4, 8, 16, 32, 64, 96, 128, 256, 512, 600, 700, 800, 900, 1024]
# [2^k for k in 0:11]

function main()

    for SimBeta in SimBetas

        plot()

        @time for N in NNBlock

            FilePathIn = PROJECT_ROOT * "/../simulations/N=$N/SimBeta=$SimBeta.txt"
            QQ = readdlm(FilePathIn, '\n', Int64; comments=true)

            QQ2 = QQ.^2
            MeanQ2 = mean(QQ2)

            @info "Blocking for different block sizes" N SimBeta MeanQ2

            ErrorsQ = fill(0.0, length(BlockSizes))
            ErrorsQ2 = fill(0.0, length(BlockSizes))

            @time for (i,k) in enumerate(BlockSizes)
                @info "Blocking Q and Q² with k=$k, number of blocks: $(round(Int64, length(QQ)/k))"
                _, ErrorsQ[i], _, ErrorsQ2[i] = BlockQ(QQ, k)
            end

            plot!(BlockSizes, ErrorsQ2, xlabel=L"k", ylabel=L"\sigma_{Q^2}",
            title=L"$\tilde\beta = %$SimBeta$", label=L"N=%$N")
        end

        savefig(PROJECT_ROOT*"/qualitative_plots/blocking/blocking_NN=$(NNBlock)_SimBeta=$SimBeta.pdf")
        println("Saved figure to /qualitative_plots/blocking/blocking_NN=$(NNBlock)_SimBeta=$SimBeta.pdf")
    end
end

@time main()
