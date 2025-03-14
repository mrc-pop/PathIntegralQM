#!/usr/bin/julia

using Printf
using DelimitedFiles
using Dates
using StatsBase
using Plots
using UnicodePlots

PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."	# Absolute path to PathIntegralQM/src
include(PROJECT_ROOT * "/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/convergence_quality.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

function WriteQVarianceSimBetaTMP(
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
				TmpEQQ = std(QMatrix[:,sb].^2, corrected=true) / sqrt(size(QMatrix,1))
				QQ[sb,1] = TmpQQ
				QQ[sb,2] = TmpEQQ
				write(DataFile, "\"$(Label)\"; $(SimBeta); $(TmpQQ); $(TmpEQQ)\n")
			end
			close(DataFile)
		end
	end
end

function PlotQVarianceSimBetaTMP(
	DirPathIn::String,
	N::Int64,
	NSweeps::Int64
)

	NSweepsString = @sprintf "%.1e" NSweeps
	DirPathOut = DirPathIn * "/convergence_plots/QVariances_SimBeta"
	mkpath(DirPathOut)
	FilePathOut = DirPathOut *
				"/N=$(N)_NSweeps=" * NSweepsString	* ".txt"

	Data = readdlm(FilePathOut, ';', comments=true)
			Labels = Data[:,1]
			SimBetas = Data[:,2]
			QQ = Data[:,3]
			eQQ = Data[:,4]

	p = plot(
		xlabel=L"$\tilde{\beta}$",
		ylabel=L"$\langle Q^2 \rangle$",
		title=L"$\chi (\tilde{\beta})$ for different local schemes ($N=%$N$)",
		legend=:topleft
	)

	LabelsUnique = unique(Labels)
	for (j,l) in enumerate(LabelsUnique)
		Indices = findall(x->x==l, Labels)
		plot!(p,
			SimBetas[Indices],  QQ[Indices],
			yerr=eQQ[Indices],
			label=l,
			markershape=:circle,
			markersize=1.5,
			linewidth=0.5,
			color=MyColors[2*j]
		)
	end

	FilePathOut = DirPathOut *
		"/N=$(N)_NSweeps=" * NSweepsString	* ".pdf"
	Plots.savefig(p, FilePathOut)

end

function main()

	DirPathIn = PROJECT_ROOT * "/../convergence"

	Waiting = true
	println("This program allows for two modes and two submodes: write and plot (w/p), deep and extended (d/e).
- Write (w): run computations on previously computed Qs, write on file;
- Plot (p): plot results;
And:
- Deep (d): do that on deep data;
- Extended (e): do that on extended data;
The four options are we, pe, wd, pd")
	print("Choose your mode: (we/pe/wd/pd) ")
	UserMode = readline()
	while Waiting
		if UserMode == "we"
			Waiting = false
			@info "Writing extended data" ConvergenceNNExtended ConvergenceSimBetasExtended

			for N in ConvergenceNNExtended
				WriteQVarianceSimBetaTMP(DirPathIn, N, NSweeps)
			end

		elseif UserMode == "pe"
			Waiting = false
			@info "Plotting extended data" ConvergenceNNExtended ConvergenceSimBetasExtended

			pgfplotsx()
			for N in ConvergenceNNExtended
				PlotQVarianceSimBetaTMP(DirPathIn, N, NSweeps)
			end

		elseif UserMode == "wd"
			Waiting = false
			@info "Writing deep data" ConvergenceNNDeep ConvergenceSimBetasDeep

			print("Mode under construction!")

		elseif UserMode == "pd"
			Waiting = false
			@info "Plotting deep data" ConvergenceNNDeep ConvergenceSimBetasDeep

			pgfplotsx()
			PlotQConvergenceFigures(DirPathIn, ConvergenceNNDeep, ConvergenceSimBetasDeep; PlotExtendedData=false, QSteps)

		else
			Waiting = true
			print("Invalid input. Please use we, pe, wd or pd to answer. Choose your mode: (we/pe/wd/pd) ")
			UserMode = readline()
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
