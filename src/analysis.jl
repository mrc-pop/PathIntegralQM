#!/usr/bin/julia

using DelimitedFiles
using Statistics
using Plots

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

const k = 1000  # block length. TODO CHANGE
const scan = false

function main()
    if length(ARGS) != 1
        @error "Wrong number of arguments. How to use this script? julia analysis.jl <SimBeta>"
        return
    end

    SimBeta = parse(Float64, ARGS[1])

    @info "Analysis parameters" SimBeta k

    if scan
        # Scan /simulations directory and find all available N
        NNScan = Int64[]
        SimDir = PROJECT_ROOT * "/../simulations/"
        for dir in readdir(SimDir)
            if !startswith(dir, "N=")
                continue # go to next dir
            end

            N = parse(Int, split(dir, "=")[2])
            FilePath = SimDir * dir * "/SimBeta=$SimBeta.txt"

            if !isfile(FilePath)
                @warn "File not found for N=$N"
                continue
            end

            push!(NNScan, N)
        end
    end

    # Arrays to store results
    Q2Means = Float64[]
    Q2Errors = Float64[]
    ττ = Float64[]

    Sizes = scan ? NNScan : NN

    for N in Sizes
        FilePath = PROJECT_ROOT * "/../simulations/N=$N/SimBeta=$SimBeta.txt"

        if !isfile(FilePath)
            @warn "File not found for N=$N"
            continue
        end

        # Read and process data
        QQ = readdlm(FilePath, '\n', Int64; comments=true)

        # Calculate quantities with blocking and store them
        _, _, Q2Mean, Q2Error = BlockQ(QQ, k)
        push!(Q2Means, Q2Mean)
        push!(Q2Errors, Q2Error)

        # Calculate tau and store it
        Q2ErrorNaive = std(QQ)  / sqrt(length(QQ))
        τ = 0.5 * ( (Q2Error / Q2ErrorNaive)^2  - 1)  # eq. (4.1.33)
        push!(ττ, τ)
    end

    ηη = SimBeta ./ NN
    ηη2 = ηη .^ 2

    RelErrorQ2 = Q2Errors ./ Q2Means

    @info "Relative error on Q²", RelErrorQ2
    @info "Tau values", ττ

    # Create plot
    p = scatter(NN, Q2Means,
            title=L"\tilde\beta = %$SimBeta",
            markersize=2,
            yerror=Q2Errors,
            xlabel=L"N",
            ylabel=L"\langle Q^2 \rangle",
            label="Data")

    # Add theoretical line (low temperature)
    if SimBeta > 1.0
        hline!([SimBeta],
            label=L"$\tilde\beta$",
            color="gray",
            linestyle=:dash)
    elseif SimBeta < 0.5
        hline!([exp(-1 / (2*SimBeta))],
            label=L"$\exp(-1/2\tilde\beta)$",
            color="gray",
            linestyle=:dash)
    end

    # Save plot
    savefig(p, PROJECT_ROOT * "/../analysis/Q2_vs_N_SimBeta=$SimBeta.pdf")
    println("Plot saved as Q2_vs_N_SimBeta=$SimBeta.pdf")

    # Using the formula Var(Q²) = 2 [τ_(Q²)]² Var(Q²)_{Naive}, calculate the autocorrelation time
    # for the observable Q². Var(Q²)_{Naive} is the one estimated with std() and without blocking.
    # cfr. Eq. (4.1.33) lecture notes
    q = scatter(NN, ττ,
            title=L"\tilde\beta = %$SimBeta",
            markersize=2,
            #yerror=Q2Errors,
            yscale=:log10,
            xlabel=L"N",
            ylabel=L"\tau",
            label="Data")

    # Save plot
    savefig(q, PROJECT_ROOT * "/../analysis/tau_vs_N_SimBeta=$SimBeta.pdf")
    println("Plot saved as tau_vs_N_SimBeta=$SimBeta.pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
