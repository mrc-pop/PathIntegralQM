#!/usr/bin/julia

"""
Plot the path given by Config.Lattice, including \"pacman effect\".
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

function PlotPathUnicode(Config::Configuration)
    Q = round(Int, CalculateQ(Config))

    PacmanUpIndices = findall(diff(Config.Lattice) .< -0.5)
    PacmanDownIndices = findall(diff(Config.Lattice) .> 0.5)

    # Add x_{N+1} = x_0 for plotting PBC explicitly
    LatticePlot = vcat(Config.Lattice, Config.Lattice[1])

    # title = "β̃ = $(round(Config.SimBeta, digits=2)), N = $(Config.N), η = $(round(Config.Eta, digits=2)), Q = $Q"

    p = lineplot(
        [1],#1:Config.N+1,
        [0],#LatticePlot,
        title = "Q = $Q",
        xlabel = "i",
        ylabel = "xᵢ",
        width = 120,
        height = 15,
        xlim = (1, Config.N+1),
        ylim = (0.0, 1.0)
    )

    for i in 1:Config.N
        if i in PacmanUpIndices
            lineplot!(p, [i, i+1], [LatticePlot[i], LatticePlot[i+1]+1], color=:blue)
            lineplot!(p, [i, i+1], [LatticePlot[i]-1, LatticePlot[i+1]], color=:blue)
        elseif i in PacmanDownIndices
            lineplot!(p, [i, i+1], [LatticePlot[i], LatticePlot[i+1]-1], color=:red)
            lineplot!(p, [i, i+1], [LatticePlot[i]+1, LatticePlot[i+1]], color=:red)
        else
            lineplot!(p, [i, i+1], [LatticePlot[i], LatticePlot[i+1]], color=:white)
        end
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
    if iEnd === nothing
        println("No Tailor plot produced, since iEnd is nothing.")
        return
    end

    i = Site
    plot!(xlims=(1,Config.N+1))
    hline!([Config.Lattice[i]], color="gray", alpha=0.3, label=nothing)

    # Find sites which should make a pacman effect
    PacmanUpIndices = findall(diff(xxNewPlot) .< -0.5)
    PacmanDownIndices = findall(diff(xxNewPlot) .> 0.5)

    # routine identical to PlotPath
    if i<iEnd
        PacmanUpIndices = findall(diff(xxNewPlot) .< -0.5)
        PacmanDownIndices = findall(diff(xxNewPlot) .> 0.5)
        for j in i:iEnd
            index = j-i+1 # varies between 1:(DeltaiEnd+1)
            if index in PacmanUpIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]+1], color="gray",
                    label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]-1, xxNewPlot[index+1]], color="gray",
                    label=nothing, linestyle=:dash)
            elseif index in PacmanDownIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]-1], color="gray",
                    label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]+1, xxNewPlot[index+1]], color="gray",
                    label=nothing, linestyle=:dash)
            else
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]], color="gray",
                    label=nothing, linestyle=:dash)
            end
        end
    else
        # Collect all points which form the cluster (from which a dashed line starts)
        ClusterPoints = vcat(i:N, 1:iEnd)
        # Plot xxNewPlot as above, since xxNewPlot[i] is the value of x
        # in the site contained in ClusterPoints[i]
        for (index,j) in enumerate(ClusterPoints)
            if index in PacmanUpIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]+1], color="gray",
                    label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]-1, xxNewPlot[index+1]], color="gray",
                    label=nothing, linestyle=:dash)
            elseif index in PacmanDownIndices
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]-1], color="gray",
                    label=nothing, linestyle=:dash)
                plot!([j, j+1], [xxNewPlot[index]+1, xxNewPlot[index+1]], color="gray",
                    label=nothing, linestyle=:dash)
            else
                plot!([j, j+1], [xxNewPlot[index], xxNewPlot[index+1]], color="gray",
                    linestyle=:dash, label=nothing)
            end
        end
    end

    # Add indicator in iEnd
    scatter!([iEnd], [mod(2*Config.Lattice[i] - Config.Lattice[iEnd],1)], color="gray",
        markersize=2, label=nothing)

    # Add dots to the start and end of cluster
    LastPoint = iEnd < Config.N ? iEnd+1 : 2
    scatter!([i, LastPoint], [Config.Lattice[i], Config.Lattice[LastPoint]], color="black",
        markersize=2, label=nothing)
end

"""
Warning function if we try to plot the cluster but iEnd is `nothing`.
"""
function PlotTailorUpdate!(
    Config::Configuration,
    Site::Int64,
    iEnd::Nothing,
    xxNewPlot::Nothing)
    println("No Tailor plot produced, since iEnd is nothing.")
    return
end

# --------------------------------- Q plots ------------------------------------

function PlotQHistogramsUnicode(
	FilePathIn::String;
	Bins=-6.5:1.0:6.5
)
	
	Data = readdlm(FilePathIn, ';', comments=true)
	SimBetas = Data[1,:]
	QMatrix = Data[2:end,:]
	
	Indices = vcat("Index", [x for x in 1:length(SimBetas)])
	Values = vcat("SimBeta", SimBetas)
	UserSelectionMatrix = hcat(Indices, Values)
	@info "Choose SimBeta to plot (enter index)" UserSelectionMatrix
	print("Choose index: (Int) ")
	UserSelection = readline()
	UserIndex = parse(Int64, UserSelection)
	SimBeta = SimBetas[UserIndex]
	
	h = fit(Histogram, QMatrix[:, UserIndex], Bins)
	h = normalize(h; mode=:pdf)
	
	unicodeplots()
	# pgfplotsx()
	
	
	xmax = 7
	xmin = -xmax
	
	p = plot(
		xlabel="Q",
		ylabel="Normalized occurrencies",
		title="Q occurrencies density",
		xlim=(xmin,xmax),
		ylim=(0,1)
	)
	
	Centers = [(Bins[i]+Bins[i+1])/2 for i in 1:length(Bins)-1]
	bar!(Centers, h.weights,
		label="Normalized distribution of Qs", 
		color=:red)
	
	xx = [x for x in xmin:1.0:xmax]

	plot!(p, xx, exp.(-xx.^2 ./ (2*SimBeta)) ./ sqrt(2*pi*SimBeta),
		label="Theoretical gaussian distribution",
		color=:blue)
	
	@info "Data from: $FilePathIn" p
	# savefig(p, "tmp.pdf")
end
