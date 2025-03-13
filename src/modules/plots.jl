#!/usr/bin/julia

using Printf

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

# ------------------------------ Q histograms ----------------------------------

function PlotQHistogramsUnicode(
	FilePathIn::String;
	Bins=-6.5:1.0:6.5
)

	Sequential, N, Scheme, NSweeps = ProcessQDataFilePath(FilePathIn)

	@info "Loading data."
	Data = readdlm(FilePathIn, ';', comments=true)
	SimBetas = Data[1,:]
	QMatrix = Data[2:end,:]
	@info "Data loaded."

	Indices = vcat("Index", [x for x in 1:length(SimBetas)], 0)
	Values = vcat("SimBeta", SimBetas, "Exit")
	UserSelectionMatrix = hcat(Indices, Values)
	Repeat = true

	while Repeat
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
			title="Q pdf; N=$N; SimBeta=$SimBeta",
			xlim=(xmin,xmax),
			ylim=(0,1)
		)

		Centers = [(Bins[i]+Bins[i+1])/2 for i in 1:length(Bins)-1]
		bar!(Centers, h.weights,
			label="Normalized distribution of Qs",
			color=:red)

		xx = [x for x in xmin:0.1:xmax]

		plot!(p, xx, exp.(-xx.^2 ./ (2*SimBeta)) ./ sqrt(2*pi*SimBeta),
			label="Theoretical gaussian distribution",
			color=:blue)

		@info "Data from: $FilePathIn" p

		print("Would you like to plot another? (y/n) ")
		UserRepeat = readline()
		Waiting = true

		while Waiting
			if UserRepeat == "y"
				Waiting = false
				Repeat = true
			elseif UserRepeat == "n"
				Waiting = false
				Repeat = false
				println("Exiting.")
			else
				Waiting = true
				print("Invalid input. Please use y or n to answer. Would you like to plot another? (y/n) ")
			end
		end
	end
end

function PlotQHistograms(
	FilePathIn::String,
	Indices::Vector{Int64};
	Bins=-6.5:1.0:6.5
)

	Sequential, N, Scheme, NSweeps = ProcessQDataFilePath(FilePathIn)

	@info "Loading data" Scheme Sequential N NSweeps
	Data = readdlm(FilePathIn, ';', comments=true)
	SimBetas = Data[1,:]
	QMatrix = Data[2:end,:]

	for Index in Indices

		SimBeta = SimBetas[Index]

		h = fit(Histogram, QMatrix[:, Index], Bins)
		h = normalize(h; mode=:pdf)

		pgfplotsx()

		xmax = 7
		xmin = -xmax

		p = plot(
			xlabel=L"$Q$",
			ylabel="Normalized occurrencies",
			title=L"$Q$ pdf ($N=%$N$, $\tilde{\beta}=%$SimBeta$)",
			xlim=(xmin,xmax),
			ylim=(0,1)
		)

		Centers = [(Bins[i]+Bins[i+1])/2 for i in 1:length(Bins)-1]
		bar!(Centers, h.weights,
			label=L"$Q$ pdf",
			color=:red)

		xx = [x for x in xmin:0.1:xmax]

		plot!(p, xx, exp.(-xx.^2 ./ (2*SimBeta)) ./ sqrt(2*pi*SimBeta),
			label="Theory",
			color=:blue)

		FilePathOut = FilePathIn[1:end-4] * "_SimBeta=$(SimBeta).pdf"
		Plots.savefig(p, FilePathOut)
	end
end

function PlotQHistogramsCompared(
	DirPathIn::String, # Up to /convergence/
	SimBeta::Float64,
	N::Int64,
	NSweeps::Int64;
	Bins=-6.5:1.0:6.5,
	BarWidth = 0.2
)

	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]

	NSweepsString = @sprintf "%.1e" NSweeps
	Histograms = Any[]
	Labels = ["MS", "MR", "HS", "HR"]

	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)
			FilePathIn = DirPathIn *
				"/" * Mode *
				"/N=$N/" *
				Scheme *
				"_NSweeps=" * NSweepsString *
				".txt"

			@info "Loading data" Scheme Mode N NSweeps SimBeta
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[1,:]
			QMatrix = Data[2:end,:]

			Index = findall(x->x==SimBeta, SimBetas)
			h = fit(Histogram, vec(QMatrix[:, Index]), Bins)
			h = normalize(h; mode=:pdf)
			push!(Histograms, [h, Labels[2*(s-1)+m]])
		end
	end

	Edges = [x for x in Bins]
	Centers = [round( (Edges[i+1]+Edges[i])/2, digits=2 ) for i in 1:length(Edges)-1]

	pgfplotsx()

	xmax = 7
	xmin = -xmax

	p = plot(
		xlabel=L"$Q$",
		ylabel="Normalized occurrencies",
		title=L"$Q$ pdf ($N=%$N$, $\tilde{\beta}=%$SimBeta$)",
		xlim=(xmin,xmax),
		ylim=(0,1),
		xticks=-6:1:6,
		minorticks=false,
		legend=:topright
	)

	for (i,k) in enumerate([1.5, 0.5, -0.5, -1.5])	# Inverted sequence
		bar!(
			Centers.-k*BarWidth, Histograms[i][1].weights,
			bar_width=BarWidth,
			color=MyColors[2*i],
			label=Histograms[i][2]
		)
	end

	xx = [x for x in xmin:0.1:xmax]

	plot!(p, xx, exp.(-xx.^2 ./ (2*SimBeta)) ./ sqrt(2*pi*SimBeta),
		label="Theory",
		color=:blue)

	DirPathOut = DirPathIn * "/convergence_plots/QHistograms_compared"
	mkpath(DirPathOut)
	FilePathOut = DirPathOut *
		"/N=$(N)" *
		"_NSweeps=" * NSweepsString *
		"_SimBeta=$(SimBeta).pdf"
	Plots.savefig(p, FilePathOut)
end

# -------------------------------- Q correlators -------------------------------

function PlotQCorrelators(
	DirPathIn::String,
	SimBeta::Float64,
	N::Int64,
	NSweeps::Int64;
	Skip=4
)

	NSweepsString = @sprintf "%.1e" NSweeps

	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]
	Labels = ["MS", "MR", "HS", "HR"]

	p = plot(
		xlabel=L"$k/n_Q$",
		ylabel=L"$C_Q(k)$",
		title=L"$Q$ correlator ($N=%$N, \tilde\beta=%$SimBeta$)",
		legend=:topright
	)

	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)
			FilePathIn = DirPathIn *
				"/" * Mode *
				"/N=$N/Q_deep/" *
				Scheme *
				"_NSweeps=" * NSweepsString *
				"_QCorrelators.txt"

			@info "Loading data" Scheme Mode N NSweeps
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[1,:]
			CMatrix = Data[2:Skip:end,:]

			Index = findall(x->x==SimBeta, SimBetas)
			CC = CMatrix[:, Index]
			kk = 1:Skip:(size(Data,1)-1)

			plot!(p, kk,  CC,
				label=Labels[2*(s-1)+m],
				markershape=:circle,
				markersize=1.5,
				linewidth=0.5,
				color=MyColors[2*(2*(s-1)+m)])
		end
	end

	DirPathOut = DirPathIn * "/convergence_plots/QCorrelators"
	mkpath(DirPathOut)
	FilePathOut = DirPathOut *
		"/N=$(N)_NSweeps=" * NSweepsString	* "_SimBeta=$SimBeta.pdf"
	Plots.savefig(p, FilePathOut)

end

# ---------------------------------- Q blocks ----------------------------------

function ParseArray(
	WeightsStr::String,
	BinsStr::String
)

	ParsedWeights = parse.(Float64, split(strip(WeightsStr, ['[', ']', ' ']), ','))

	BinsStr = strip(BinsStr, ' ')
	Delims = findall(':', BinsStr)
	Start = parse(Float64, BinsStr[1:Delims[1]-1])
	Step = parse(Float64, BinsStr[Delims[1]+1:Delims[2]-1])
	Stop = parse(Float64, BinsStr[Delims[2]+1:end])
	ParsedBins = Start:Step:Stop
	return ParsedWeights, ParsedBins
end

function PlotQBlocksHistogramsCompared(
	DirPathIn::String, # Up to /convergence/
	SimBeta::Float64,
	N::Int64,
	NSweeps::Int64;
	BarWidth = 0.2
)

	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]

	NSweepsString = @sprintf "%.1e" NSweeps
	Histograms = Any[]
	Labels = ["MS", "MR", "HS", "HR"]
	Spacings = [0.5, -0.5]

	pgfplotsx()

	xmax = 31
	xmin = -1

	for (s,Scheme) in enumerate(Schemes)

		p = plot(
			xlabel="Block length",
			ylabel="Occurrencies",
			title=L"Same-$Q$ block lengths (%$Scheme, $N=%$N$, $\tilde{\beta}=%$SimBeta$)",
			xlim=(xmin,xmax),
			ylim=(0,1),
			xticks=0:10:100,
			minorticks=false,
			legend=:topright,
			size=(700,300)
		)

		for (m,Mode) in enumerate(Modes)
			FilePathIn = DirPathIn *
				"/" * Mode *
				"/N=$N/Q_deep/" *
				Scheme *
				"_NSweeps=" * NSweepsString *
				"_QHistograms.txt"

			@info "Loading data" Scheme Mode N NSweeps SimBeta
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[:,1]
			Index = findall(x->x==SimBeta, SimBetas)[1]

			Weights, Bins = ParseArray(String(Data[Index, 2]), String(Data[Index, 3]))

			h = fit(Histogram, Weights, Bins)
			h = normalize(h; mode=:pdf)
			push!(Histograms, [h, Labels[2*(s-1)+m]])

			Edges = [x for x in Bins]
			Centers = [round( (Edges[i+1]+Edges[i])/2, digits=2 ) for i in 1:length(Edges)-1]

			Index = 2*(s-1)+m	# Overwrite variable

			bar!(
				Centers.-Spacings[m]*BarWidth, h.weights,
				bar_width=BarWidth,
				color=MyColors[2*Index],
				label=Labels[Index]
			)

		end

		DirPathOut = DirPathIn * "/convergence_plots/QBlocksHistograms_compared"
		mkpath(DirPathOut)
		FilePathOut = DirPathOut *
			"/N=$(N)" *
			"$(Scheme)_NSweeps=" * NSweepsString *
			"_SimBeta=$(SimBeta).pdf"
		Plots.savefig(p, FilePathOut)
	end
end


# --------------------------------- Q variance ---------------------------------

function PlotQVarianceSimBeta(
	DirPathIn::String,
	N::Int64,
	NSweeps::Int64
)

	NSweepsString = @sprintf "%.1e" NSweeps
	DirPathOut = DirPathIn * "/convergence_plots/QVariances_SimBeta"
	mkpath(DirPathOut)
	FilePathOut = DirPathOut *
				"/N=$(N)_NSweeps=" * NSweepsString	* ".txt"

	DataFile = open(FilePathOut, "w")
	write(DataFile, "# Label; SimBeta; <Q^2>; eQ^2\n")
	close(DataFile)

	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]
	Labels = ["MS", "MR", "HS", "HR"]

	p = plot(
		xlabel=L"$\tilde{\beta}$",
		ylabel=L"$\langle Q^2 \rangle$",
		title=L"$\chi (\tilde{\beta})$ for different local schemes ($N=%$N$)",
		legend=:topleft
	)

	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)
			Label = Labels[2*(s-1)+m]

			FilePathIn = DirPathIn *
				"/" * Mode *
				"/N=$N/" *
				Scheme *
				"_NSweeps=" * NSweepsString *
				".txt"

			@info "Loading data" Scheme Mode N NSweeps
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[1,:]
			QMatrix = Data[2:end,:]

			DataFile = open(FilePathOut, "a")
			QQ = zeros(length(SimBetas),2)
			for (sb, SimBeta) in enumerate(SimBetas)
				TmpQQ = mean(QMatrix[:,sb].^2)
				TmpEQQ = std(QMatrix[:,sb].^2)
				QQ[sb,1] = TmpQQ
				QQ[sb,2] = TmpEQQ
				write(DataFile, "$(Label); $(SimBeta); $(TmpQQ); $(TmpEQQ)\n")
			end
			close(DataFile)

			plot!(p, SimBetas,  QQ[:,1],
				label=Label,
				markershape=:circle,
				markersize=1.5,
				linewidth=0.5,
				color=MyColors[2*(2*(s-1)+m)])
		end
	end

	FilePathOut = DirPathOut *
		"/N=$(N)_NSweeps=" * NSweepsString	* ".pdf"
	Plots.savefig(p, FilePathOut)

end

# ------------------------- Convergence plots handler --------------------------

function PlotQConvergenceFigures(DirPathIn, CNN, CSimBetas; PlotExtendedData=false)

	for N in CNN

		println("Plotting N=$N data")
		if !PlotExtendedData
			for SimBeta in CSimBetas
				PlotQHistogramsCompared(DirPathIn, SimBeta, N, NSweeps)
				PlotQCorrelators(DirPathIn, SimBeta, N, NSweeps)
			end
		elseif PlotExtendedData
			PlotQVarianceSimBeta(DirPathIn, N, NSweeps)
		end
		println()

	end
end
