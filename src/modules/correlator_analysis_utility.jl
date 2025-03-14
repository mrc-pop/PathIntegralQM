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

function WriteQDataDiffs(
	DirPathIn::String,	# up to /convergence/
	N::Int64,
	NSweeps::Int64
)
	
	Now = now()	# Comparable timestamps

	# Prepare files headers: first line contains SimBetas
	SimBetasHeader = "# First line: SimBetas\n"
	QDiffsHeader = "# Following lines: QDiffs (divided for each SimBeta) [calculated $Now]\n"

	NSweepsString = @sprintf "%.1e" NSweeps
	DirPathOut = DirPathIn * "/convergence_plots/QDiffs"
	mkpath(DirPathOut)
	
	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]
	
	DiffsDict = Dict()
	
	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)

			FilePathIn = DirPathIn *
				"/" * Mode *
				"/N=$N/" *
				Scheme *
				"_NSweeps=" * NSweepsString *
				".txt"

			@info "Loading data" Scheme Mode N NSweeps
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[1,:]
			
			SimBetasString = ""
			for SimBeta in SimBetas
				SimBetasString *= "$SimBeta; "
			end
			SimBetasChars = collect(SimBetasString)
			popat!(SimBetasChars, length(SimBetasString)-1)
			SimBetasString = String(SimBetasChars) * "\n"
			
			QMatrix = Data[2:end,:]		
			DiffMatrix = zeros(size(QMatrix,1)-1, size(QMatrix,2))
			
			for (sb,SimBeta) in enumerate(SimBetas)
				QQ = QMatrix[:,sb]
				DiffMatrix[:,sb] .= diff(QQ)
			end
			
			DiffsDict[Scheme * "-" * Mode] = DiffMatrix
			
			FilePathOut = DirPathOut *
				"/" * Scheme * "_" * Mode * 
				"N=$(N)_NSweeps=" * NSweepsString * ".txt"

			DataFile = open(FilePathOut, "w")
			write(DataFile, SimBetasHeader)
			write(DataFile, SimBetasString)
			write(DataFile, QDiffsHeader)
			close(DataFile)
			
			open(FilePathOut, "a") do io
				writedlm(io, DiffMatrix, "; ")
			end
		end
	end
	
	return DiffsDict
	
end

function ReadQDataDiffs(
	DirPathIn::String,	# up to /convergence/
	N::Int64,
	NSweeps::Int64
)
	
	NSweepsString = @sprintf "%.1e" NSweeps
	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]
	
	DiffsDict = Dict()
	
	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)

			FilePathIn = DirPathIn *
				"/convergence_plots/QDiffs" *
				"/N=$(N)_NSweeps=" * NSweepsString	* ".txt"

			@info "Loading data" Scheme Mode N NSweeps
			Data = readdlm(FilePathIn, ';', comments=true)
			SimBetas = Data[1,:]
			DiffMatrix = Data[2:end,:]
			DiffsDict[Scheme * "-" * Mode] = DiffMatrix
		end
	end
	
	return DiffsDict
	
end

function PlotQDataDiffs(
	DirPathIn::String,	# up to /convergence/
	N::Int64,
	NSweeps::Int64;
	kMax = 100
)
	
	NSweepsString = @sprintf "%.1e" NSweeps
	Schemes = ["Metropolis", "Heatbath"]
	Modes = ["sequential", "random"]
	Labels = ["MS", "MR", "HS", "HR"]
	
	DiffsDict = Dict()
	
	p = plot(
		xlabel=L"$k$",
		ylabel=L"$D(k)$",
		title=L"$D(k)$ for different local schemes ($N=%$N$)",
		legend=:topleft
	)
	
	DiffDict = ReadQDataDiffs(DirPathIn, N, NSweeps)
	
	for (s,Scheme) in enumerate(Schemes)
		for (m,Mode) in enumerate(Modes)
			
			DiffMatrix = DiffDict[Scheme * "-" * Mode]
			
			kk = 1:1:kMax
			DD = zeros(kMax)
			L = size(DiffMatrix,1)
			
			for k in kk
			
				TmpD = 0
				for j in 1:L-k
					TmpD += abs(DiffMatrix[j]+DiffMatrix[j+k])
				end
				TmpD /= L-k
				
				DD[k] = TmpD
			end
			
			plot!(p,
				kk, DD,
				label=Labels[2*(s-1)+m],
				markershape=:circle,
				markersize=1.5,
				linewidth=0.5,
				color=MyColors[2*(2*(s-1)+m)]
			)
			
		end
	end
	
	gui()
	
	return DiffsDict
	
end

function main()
	
	DirPathIn = PROJECT_ROOT * "/../convergence"
	DiffsDict = ReadQDataDiffs(DirPathIn, 50, NSweeps)
	UserMode = "Heatbath-sequential"
	
	sb = 2
	DD = DiffsDict[UserMode]
	kk = 5:5:50

	for k in kk
		Counter = 0
		for j in 1:size(DD,1)-k
			if DD[j,sb]+DD[j+k,sb]==0
				Counter += 1
			end
		end
		println("Competing diffs ratio at k=$k: $(Counter/size(DD,1))")
	end
	
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
