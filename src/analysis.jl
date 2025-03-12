#!/usr/bin/julia

using DelimitedFiles
using Statistics
using Plots

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/../src/modules/processing.jl")

# const LengthsBlockSizes = Dict( # see txt file
#     # N => k
#     75 => 10,
#     100 => 25,
#     125 => 50,
#     150 => 100,
#     175 => 100,
#     200 => 100,
#     300 => 250,
#     400 => 1000
# )

const LengthsBlockSizes = Dict( # see txt file
    # N => k
    75 => 50,
    100 => 50,
    125 => 50,
    150 => 50,
    175 => 50,
    200 => 50,
    300 => 100,
    400 => 100
)

# const LengthsBlockSizes = Dict(
#     # N => k
#     50 => 1000,
#     75 => 1000,
#     100 => 1000,
#     125 => 1000,
#     150 => 1000,
#     175 => 1000,
#     200 => 1000,
#     300 => 1000,
#     400 => 10000
# )


const Scan = true          # Scan all available N # TODO remove
const MakeHist = false      # plot Histogram

function main()
    if length(ARGS) != 1
        @error "Wrong number of arguments. How to use this script? julia analysis.jl <SimBeta>"
        return
    end

    # Define kk as the values of LengthsBlockSizes
    Sizes = collect(keys(LengthsBlockSizes))
    kk = collect(values(LengthsBlockSizes))

    SimBeta = parse(Float64, ARGS[1])

    @info "Analysis parameters" SimBeta LengthsBlockSizes Scan

    if Scan
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
    QMeans = Float64[]
    QErrors = Float64[]
    Q2Means = Float64[]
    Q2Errors = Float64[]
    ττ = Float64[]

    Sizes = Scan ? NNScan : Sizes

    for (SizeIndex, N) in enumerate(Sizes)

        @info "Performing blocking on N=$N..."

        FilePath = PROJECT_ROOT * "/../simulations/N=$N/SimBeta=$SimBeta.txt"

        if !isfile(FilePath)
            @warn "File not found for N=$N"
            continue
        end

        # Read and process data
        QQ = readdlm(FilePath, '\n', Int64; comments=true)

        if MakeHist
            println("Unique Q values for N=$N, $(unique(QQ))")
            # Histogram plot
            mkpath(PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/hist/", )
            bins = collect(minimum(QQ)-0.5:maximum(QQ)+0.5)
            hh = plot(QQ, seriestype=:barhist, bins=bins, xlims=(-12, 12), ylims=(0.0,0.3),
                normalize=:pdf, xlabel=L"Q", ylabel="rel. counts", label=:none,
                title=L"\tilde \beta = %$SimBeta, N = %$N")
            savefig(hh, PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/hist/hist_N=$N.pdf")
            println("Plot saved as /SimBeta=$SimBeta/hist/hist_N=$N.pdf")
        end
        # Calculate quantities with blocking and store them
        QMean, QError, Q2Mean, Q2Error = BlockQ(QQ, kk[SizeIndex])
        push!(QMeans, QMean)
        push!(QErrors, QError)
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

    println()

    @info "Relative error on Q²", RelErrorQ2
    @info "Tau values", ττ

    println()

    mkpath(PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/", )

    # Plot name
    if TailorSteps[1] == TailorSteps[end]
        TS = "TailorSteps=$(TailorSteps[1])"
    else
        TS = "TailorSteps=$((TailorSteps./NN)[end])NN"
    end

    if NR !== 1
        PT = "NR_$NR"
    end

    # Create plot
    p = scatter(Sizes, Q2Means,
            title=L"\tilde\beta = %$SimBeta",
            markersize=2,
            yerror=Q2Errors,
            xlabel=L"N",
            ylabel=L"\langle Q^2 \rangle",
            label=:none
        )

    # Add theoretical line
    if SimBeta > 1.0
        hline!([SimBeta],
            label=:none,#L"$\tilde\beta$",
            color="gray",
            linestyle=:dash)
    elseif SimBeta < 0.5
        hline!([exp(-1 / (2*SimBeta))],
            label=:none,#L"$\exp(-1/2\tilde\beta)$",
            color="gray",
            linestyle=:dash)
    end

    # Save plot
    savefig(p, PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/Q2_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")
    println("Plot saved as SimBeta=$SimBeta/Q2_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")

    # Create plot of Q vs NN
    p2 = scatter(Sizes, QMeans,
            title=L"\tilde\beta = %$SimBeta",
            markersize=2,
            yerror=QErrors,
            xlabel=L"N",
            ylabel=L"\langle Q \rangle",
            label=:none)
    hline!([0.0],
            label=:none,
            color="gray",
            linestyle=:dash)
    # Save plot
    savefig(p2, PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/Q_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")
    println("Plot saved as SimBeta=$SimBeta/Q_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")


    # Using the formula Var(Q²) = (1 + 2 [τ_(Q²)]) Var(Q²)_{Naive}, calculate the autocorrelation time
    # for the observable Q². Var(Q²)_{Naive} is the one estimated with std() and without blocking.
    # cfr. Eq. (4.1.33) lecture notes
    q = scatter(Sizes, ττ,
            title=L"\tilde\beta = %$SimBeta",
            markersize=2,
            #yerror=Q2Errors,
            #yscale=:log10,
            xlabel=L"N",
            ylabel=L"\tau",
            label=:none)

    savefig(q, PROJECT_ROOT * "/../analysis/SimBeta=$SimBeta/tau_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")
    println("Plot saved as SimBeta=$SimBeta/tau_vs_N_SimBeta=$(SimBeta)_"*TS*"_NR=$NR.pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
