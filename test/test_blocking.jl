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
const BlockSizes = [1,5,10,15,20,40,60,80,100,200,300,400]

function main()

    for SimBeta in SimBetas

        @time for N in NN

            QQ2 = QQ.^2
            MeanQ2 = mean(QQ2)

            @info "Blocking for different block sizes" N SimBeta MeanQ2

            FilePathIn = PROJECT_ROOT * "/../simulations/N=$N/SimBeta=$SimBeta.txt"
            QQ = readdlm(FilePathIn, '\n', Int64; comments=true)

            ErrorsQ = fill(0.0, length(BlockSizes))
            ErrorsQ2 = fill(0.0, length(BlockSizes))

            for (i,k) in enumerate(BlockSizes)
                @info "Blocking Q and Q² with k=$k, number of blocks: $(round(Int64, length(QQ)/k))"
                _, ErrorsQ[i], _, ErrorsQ2[i] = BlockQ(QQ, k)
            end

            plot!(BlockSizes, ErrorsQ2, xlabel=L"$k$", ylabel=L"$\sigma_{Q^2}$",
                title=L"$\tilde \beta = %$(round(SimBeta, digits=2)),
                    N_\mathrm{sweeps} = 10^8$", label=L"$N=%$N$")
        end

        savefig(PROJECT_ROOT*"/qualitative_plots/blocking_NN=$(NN)_SimBeta=$SimBeta.pdf")
    end
end

@time main()
