#!/usr/bin/julia

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/montecarlo.jl")
include(PROJECT_ROOT * "/../src/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/../src/modules/plots.jl")

using Statistics

# const N = 100
const SimBeta = 2.0
const Δ = 0.5
const NSweeps = 5000
const NSweepsTherm = 10000

function main1()

    println("""
    Test of Metropolis updates and plot of the resulting path; tailor update and plot
    of the resulting path, with an indication of whether it has been accepted.
    """)

    Config = SetLattice(SimBeta, N)

    for i in 1:(10*N)
        Site = mod1(i,N)
        MetropolisUpdate!(Config, Site; Δ)
    end

    println("Before tailor, Q=$(CalculateQ(Config))")

    p = PlotPath(Config)

    savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$N.pdf")

    plot(p)

    Site = 10
    Found, Acc, iEnd, xxNewPlot = TailorUpdate!(Config, Site, 0.02)

    println("iEnd=$iEnd. Found? $(Bool(Found)). Accepted? $(Bool(Acc))")

    PlotTailorUpdate!(Config, Site, iEnd, xxNewPlot)

    println("After tailor, Q=$(CalculateQ(Config))")

    if Bool(Found) == 1
        savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$(N)_tailor.pdf")
        PlotPath(Config)
        savefig(PROJECT_ROOT * "/qualitative_plots/path_N=$(N)_after.pdf")
    end

end

"""
Plot Q vs MC time, measuring once every QSweepStep lattice sweeps, for the chosen sizes.
"""
function PlotQvsTime(NN::Vector{Int64}, QSweepStep::Int64; FilePathOut="")

    @info "Parameters for time plot" NN NSweeps

    QQ_Matrix = Array{Int64}[]

    for (n, N) in enumerate(NN)

        # Initialize lattice with size 10
        QStep = QSweepStep * N
        Config = SetLattice(SimBeta, N)

        # Initialize 0 matrix
        NumberOfMeasurements = round(Int,(NSweeps*N)/QStep)
        MeasureCount = 0
        QQ = fill(0, NumberOfMeasurements)

        for j in 1:(NSweeps*N)

            Site = mod1(j,N)

            MetropolisUpdate!(Config, Site; Δ)

            # Measure Q
            if mod(j, QStep) == 0
                MeasureCount += 1
                QQ[MeasureCount] = CalculateQ(Config)
            end
        end

        push!(QQ_Matrix, QQ)
    end

    z = 0.435
    plot(layout=grid(length(QQ_Matrix), 1, heights=[(1-z)/2, (1-z)/2, z]), ylims=(-4,4), grid=:none, minorticks = false)

    for (q, QQ) in enumerate(QQ_Matrix)
        if q < length(QQ_Matrix)
            plot!(QQ, label=:none, bottom_margin = -12*Plots.mm, ylabel=L"Q", subplot=q, xformatter=_->"")
            annotate!((0.02,0.8), text(L"$N=%$(NN[q])$", :left, 10), subplot=q)

        else
            plot!(1:length(QQ), QQ, subplot=q, label=:none, xlabel="Sweep", ylabel=L"Q")
            annotate!((0.015,0.8), text(L"$N=%$(NN[q])$", :left, 10), subplot=q)
        end
    end

    if FilePathOut == ""
        gui()
    else
        savefig(PROJECT_ROOT * FilePathOut)
        println("Saved time plot at " * FilePathOut)
    end
end

NN = [250, 275, 300]

PlotQvsTime(NN, 1; FilePathOut="/qualitative_plots/time/time_NN=$NN.pdf")
