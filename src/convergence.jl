#!/usr/bin/julia

using Printf
using DelimitedFiles
using Dates
using StatsBase
using Plots
using UnicodePlots

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
include(PROJECT_ROOT * "/modules/convergence_quality.jl")
include(PROJECT_ROOT * "/modules/plots.jl")

"""
Target: explore simulations quality using pure local algorithms (Metropolis and
Heatbath, in order to select the most efficient combination). Four possible
configurations:
- Sequential Metropolis
- Random selection Metropolis
- Sequential Heatbath
- Random selection Heatbath
"""
function main()

	DeepAnalysis = false
	Waiting = true
	println("This program allows for two modes: extended and deep (e/d).
- Extended (e): many NNs, many SimBetas, does not compute correlators and other deep observables;
- Deep (d): few NNs, few SimBetas, computes correlators and other deep observables.")
	print("Choose your mode: (e/d) ")
	UserMode = readline()
	while Waiting
		if UserMode == "e"
			Waiting = false
			@info "Starting extended simulation" ConvergenceNNExtended ConvergenceSimBetasExtended
			global CNN = ConvergenceNNExtended
			global CSimBetas = ConvergenceSimBetasExtended

		elseif UserMode == "d"
			Waiting = false
			@info "Starting in-depth simulation" ConvergenceNNDeep ConvergenceSimBetasDeep
			DeepAnalysis = true
			global CNN = ConvergenceNNDeep
			global CSimBetas = ConvergenceSimBetasDeep

		else
			Waiting = true
			print("Invalid input. Please use e or d to answer. Choose your mode: (e/d) ")
			UserMode = readline()
		end
	end

	Now = now()	# Comparable timestamps

	# Prepare files headers: first line contains SimBetas
	SimBetasHeader = "# First line: SimBetas\n"
	QHeader = "# Following lines: Qs (divided for each SimBeta) [calculated $Now]\n"
	QCorrelatorsHeader = "# Following lines: Q correlator at increasing distance (divided for each SimBeta) [calculated $Now]\n"
	QHistogramsHeader = "# SimBeta; Q weights; Q edges; [calculated $Now]\n"

	SimBetasString = ""
	for SimBeta in CSimBetas
		SimBetasString *= "$SimBeta; "
	end
	SimBetasChars = collect(SimBetasString)
	popat!(SimBetasChars, length(SimBetasString)-1)
	SimBetasString = String(SimBetasChars) * "\n"

	AdditionalHeaders = Dict([
		("Q", QHeader),
		("QCorrelators", QCorrelatorsHeader),
		("QHistograms", QHistogramsHeader),
		("SimBetas", SimBetasHeader),
		("SimBetasValues", SimBetasString)
	])

	NSweepsString = @sprintf "%.1e" NSweeps
	RunTimes = zeros(Float64, length(CNN), 3)

	for Seq in [true, false]

		if Seq
			DirPathOut = PROJECT_ROOT * "/../convergence/sequential/"
		elseif !Seq
			DirPathOut = PROJECT_ROOT * "/../convergence/random/"
		end

		TmpDirPathOut = DirPathOut

		for (n,N) in enumerate(CNN)

			QStep = QSteps[n]

			if Seq
				DirPathOut = TmpDirPathOut * "/N=$(N)/"
			elseif !Seq
				DirPathOut = TmpDirPathOut * "/N=$(N)/"
			end

			mkpath(DirPathOut)

			println()
			@info "Starting simulations for N=$N..."
			Results = RunConvergenceSimulations(N, NSweepsTherm, NSweeps, CSimBetas, QStep; Sequential=Seq)
			MetropolisQMatrix, MetropolisElapsedTime, HeatbathQMatrix, HeatbathElapsedTime = Results
			QMatrices = Dict([("Metropolis-Seq=$Seq", MetropolisQMatrix),
							  ("Metropolis-Seq=$Seq-Time", MetropolisElapsedTime),
							  ("Heatbath-Seq=$Seq", HeatbathQMatrix),
							  ("Heatbath-Seq=$Seq-Time", MetropolisElapsedTime)])

			RunTimes[n,:] = [N, MetropolisElapsedTime, HeatbathElapsedTime]
			@info "Total Metropolis runtime for size $N" MetropolisElapsedTime
			@info "Total Heatbath runtime for size $N" HeatbathElapsedTime

			for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])

				ElapsedTime = QMatrices[Scheme * "-Seq=$Seq-Time"]
				FilePathOut = DirPathOut * "/$(Scheme)_NSweeps=" * NSweepsString * ".txt"
				GeneralHeader = "# " * Scheme *
								", Sequential=$(Seq)," *
								" NSweeps=" * NSweepsString *
								", Elapsed time=$ElapsedTime\n"

				DataFile = open(FilePathOut, "w")
				write(DataFile, GeneralHeader)
				write(DataFile, AdditionalHeaders["SimBetas"])
				write(DataFile, AdditionalHeaders["SimBetasValues"])
				write(DataFile, AdditionalHeaders["Q"])
				close(DataFile)

				open(FilePathOut, "a") do io
					M = QMatrices[Scheme * "-Seq=$Seq"]
				    writedlm(io, M, "; ")
				end
			end

			# -------------------------- In-depth analysis ---------------------

			if DeepAnalysis

				RunDeepAnalysis(
					DirPathOut,
					NSweepsString,
					MetropolisQMatrix,		# InputMetroData
					HeatbathQMatrix,		# InputHeatData
					CSimBetas,
					Seq,
					AdditionalHeaders;
					kMax=50, # maximum number of k points
                    Skip=1
				)

				println("In-depth analysis completed.")

			end

		end

		# ------------------------------ Runtimes plot -------------------------

		unicodeplots()
		p = lineplot(
			xlabel = "N",
			ylabel = "Time [s]",
			xlim = (RunTimes[1,1], RunTimes[end,1]),
			ylim = (0, maximum(vcat(RunTimes[:,2], RunTimes[:,3]))),
			title = "Runtimes"
		)
		scatterplot!(p, RunTimes[:,1], RunTimes[:,2], color=:red, name="Metropolis", marker=:circle)
		scatterplot!(p, RunTimes[:,1], RunTimes[:,3], color=:blue, name="Heatbath", marker=:circle)
		@info "Size-wise runtimes for the following setup" CSimBetas NSweeps Seq p

		if DeepAnalysis
			FilePathOut = DirPathOut * "/../RunTimes_NSweeps=" * NSweepsString * ".txt"
			GeneralHeader = "# N; Metropolis-Seq=$(Seq); Heatbath-Seq=$(Seq)\n"

			DataFile = open(FilePathOut, "w")
			write(DataFile, GeneralHeader)
			close(DataFile)
			open(FilePathOut, "a") do io
				writedlm(io, RunTimes, "; ")
			end
		end

	end

	# ------------------------------ Detailed plot -----------------------------

	DirPathIn = PROJECT_ROOT * "/../convergence"
	pgfplotsx()
	PlotQConvergenceFigures(DirPathIn, CNN, CSimBetas; PlotExtendedData=!DeepAnalysis, QSteps)

end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
