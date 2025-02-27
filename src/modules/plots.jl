#!/usr/bin/julia

"""
Plot the path given by Config.Lattice, including "pacman effect".
"""
function PlotPath(Config::Configuration; location="nw")

    Q = round(Int, CalculateQ(Config))

    p = plot(
        xlims=(1,Config.N+1),
        ylims=(0.0,1.0),
        xlabel=L"i",
        ylabel=L"x_i",
        title=L"""
            $\tilde \beta = %$(round(Config.SimBeta, digits=2)),
            N = %$(Config.N), \eta = %$(round(Config.Eta, digits=2))$
            """
    )

    PacmanUpIndices = findall(diff(Config.Lattice) .< -0.5)
    PacmanDownIndices = findall(diff(Config.Lattice) .> 0.5)

    # Add x_{N+1} = x_0 for plotting PBC explicitly
    LatticePlot = vcat(Config.Lattice, Config.Lattice[1])

    for i in 1:Config.N
        if i in PacmanUpIndices
            plot!([i, i+1], [LatticePlot[i], LatticePlot[i+1]+1], color=MyColors[4], label=nothing)
            plot!([i, i+1], [LatticePlot[i]-1, LatticePlot[i+1]], color=MyColors[4], label=nothing)
        elseif i in PacmanDownIndices
            plot!([i, i+1], [LatticePlot[i], LatticePlot[i+1]-1], color=MyColors[1], label=nothing)
            plot!([i, i+1], [LatticePlot[i]+1, LatticePlot[i+1]], color=MyColors[1], label=nothing)
        else
            plot!([i, i+1], [LatticePlot[i], LatticePlot[i+1]], color="black", label=nothing)
        end
    end

    if location=="nw"
        annotate!((0.02,0.92), text(L"$Q=%$Q$", :left, 10))
    elseif location=="ne"
        annotate!((0.98,0.92), text(L"$Q=%$Q$", :right, 10))
    end

    return p
end

"""
Plot the tailor update as a dashed line with highlighted endpoints. Write over current plot.
xxNewPlot is the output of the TailorUpdate! function.
"""
function PlotTailorUpdate!(
    Config::Configuration,
    Site::Int64,
    iEnd::Int64,
    xxNewPlot::Array{Float64}
)
    i = Site
    plot!(xlims=(1,Config.N+1))
    hline!([Config.Lattice[i]], color="gray", alpha=0.3, label=nothing)

    PacmanUpIndices = findall(diff(xxNewPlot) .< -0.5)
    PacmanDownIndices = findall(diff(xxNewPlot) .> 0.5)

    # routine identical to PlotPath
    if i<iEnd
        for j in i:iEnd
            index = j-i+1 # varies between 1:(DeltaiEnd+1)
            if index in PacmanUpIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]+1], color="black", label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]-1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
            elseif index in PacmanDownIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]-1], color="black", label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]+1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
            else
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]], color="black", linestyle=:dash, label=nothing)
            end
        end
    else
        for j in vcat(i:N, 1:iEnd)
            # TODO FIX!!
            index = j-i+1 # varies between
            if index in PacmanUpIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]+1], color="black", label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]-1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
            elseif index in PacmanDownIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]-1], color="black", label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]+1, xxNewPlot[index+1]], color="black", label=nothing, linestyle=:dash)
            else
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]], color="black", linestyle=:dash, label=nothing)
            end
        end
    end

    scatter!([i, iEnd+1], [Config.Lattice[i], Config.Lattice[iEnd+1]], color="black",
        markersize=2, label=nothing)
end
