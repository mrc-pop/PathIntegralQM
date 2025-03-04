#!/usr/bin/julia

using Printf
using DelimitedFiles
using Dates
using StatsBase

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
	print("Should I perform a deep analysis on Q samples? (y/n) ")
	UserInput = readline()
	while Waiting
		if UserInput == "y"
			Waiting = false
			@info "Deep analysis accepted"
			DeepAnalysis = true
		elseif UserInput == "n"
			Waiting = false
			@info "Deep analysis rejected"
		else
			Waiting = true
			print("\nInvalid input. Please use y or n to answer. Should I perform a deep analysis on Q samples? (y/n) ")
		end
	end
	
	Now = now()	# Comparable timestamps

	# Prepare files headers: first line contains SimBetas
	SimBetasHeader = "# First line: SimBetas\n"
	QHeader = "# Following lines: Qs (divided for each SimBeta) [calculated $Now]\n"
	QCorrelatorsHeader = "# Following lines: Q correlator at increasing distance (divided for each SimBeta) [calculated $Now]\n"
	QHistogramsHeader = "# SimBeta; Q weights; Q edges; [calculated $Now]\n"
	
	SimBetasString = ""
	for SimBeta in SimBetas
		SimBetasString *= "$SimBeta; "
	end
	SimBetasChars = collect(SimBetasString)
	popat!(SimBetasChars, length(SimBetasString)-1)
	SimBetasString = String(SimBetasChars) * "\n"

	NSweepsString = @sprintf "%.1e" NSweeps

	for (n,N) in enumerate(NN)
	
		QStep = QSteps[n]
	
		if Sequential
			DirPathOut = PROJECT_ROOT * "/../convergence/sequential/N=$(N)/"
		elseif !Sequential
			DirPathOut = PROJECT_ROOT * "/../convergence/random/N=$(N)/"
		end
		mkpath(DirPathOut)

		println()
		@info "Starting simulations for N=$N..."
		Results = RunConvergenceSimulations(N, NSweepsTherm, NSweeps, SimBetas, QStep; Sequential)
		MetropolisQMatrix, MetropolisElapsedTime, HeatbathQMatrix, HeatbathElapsedTime = Results
		QMatrices = Dict([("Metropolis-Seq=$Sequential", MetropolisQMatrix),
						  ("Metropolis-Seq=$Sequential-Time", MetropolisElapsedTime), 
						  ("Heatbath-Seq=$Sequential", HeatbathQMatrix),
						  ("Heatbath-Seq=$Sequential-Time", HeatbathElapsedTime)])

		for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])

			ElapsedTime = QMatrices[Scheme * "-Seq=$Sequential-Time"]
			FilePathOut = DirPathOut * "/$(Scheme)_NSweeps=$(NSweepsString).txt"
			GeneralHeader = "# " * Scheme *
							", Sequential=$(Sequential)," *
							" NSweeps=" * NSweepsString *
							", Elapsed time=$ElapsedTime\n"
			
			DataFile = open(FilePathOut, "w")
			write(DataFile, GeneralHeader)
			write(DataFile, SimBetasHeader)
			write(DataFile, SimBetasString)
			write(DataFile, QHeader)
			close(DataFile)
			
			open(FilePathOut, "a") do io
				M = QMatrices[Scheme * "-Seq=$Sequential"]
		        writedlm(io, M, "; ")
		    end
		end
		
		# -------------------------- In-depth analysis -------------------------
		
		# TODO Factorize: write an algorithm that can also read data files not
		# TODO generate them from scratch each time.
		
		if DeepAnalysis
			kMax = 100
			Bins = -5.5:1.0:5.5 # -5, -4, ..., 4, 5
			
			MetropolisQCorrelators = zeros(kMax+1,length(SimBetas))
			HeatbathQCorrelators = zeros(kMax+1,length(SimBetas))
			QCorrelators = Dict()
			
			MetropolisQBlockLengths = Any[]
			HeatbathQBlockLengths = Any[]
			
			MetropolisHistogram = Any[]	# TODO SISTEMA (anche linee 134-135)
			HeatbathHistogram = Any[]	# TODO SISTEMA (anche linee 134-135)
			QHistograms = Dict()
			
			for (sb,SimBeta) in enumerate(SimBetas)
				for k in 1:kMax+1
					MetropolisQCorrelators[k,sb] = GetQCorrelator(k, MetropolisQMatrix[:,sb])
					HeatbathQCorrelators[k,sb] = GetQCorrelator(k, HeatbathQMatrix[:,sb])
				end
									
				MetroQBL = GetQBlockLengths(MetropolisQMatrix[:,sb])
				push!(MetropolisQBlockLengths, MetroQBL)
				
				HeatQBL = GetQBlockLengths(HeatbathQMatrix[:,sb])
				push!(HeatbathQBlockLengths, HeatQBL)
			
				hM = fit(Histogram, MetroQBL, Bins)
				push!(MetropolisHistogram, [SimBeta, hM.weights, hM.edges])
				
				hH = fit(Histogram, HeatQBL, Bins)
				push!(HeatbathHistogram, [SimBeta, hH.weights, hH.edges])
			
			end
			
			QCorrelators = Dict([("Metropolis-Seq=$Sequential", MetropolisQCorrelators), 
								("Heatbath-Seq=$Sequential", HeatbathQCorrelators)])
			
			QHistograms = Dict([("Metropolis-Seq=$Sequential", MetropolisHistogram), 
								 ("Heatbath-Seq=$Sequential", HeatbathHistogram)])
			
			for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])
				mkpath(DirPathOut * "/Q_deep/")
				
				# Write correlators on file
				FilePathOut = DirPathOut * "/Q_deep/$(Scheme)_NSweeps=$(NSweepsString)_QCorrelators.txt"
				GeneralHeader = "# " * Scheme * ", Sequential=$(Sequential), NSweeps=" * NSweepsString * "\n"
				
				DataFile = open(FilePathOut, "w")
				write(DataFile, GeneralHeader)
				write(DataFile, SimBetasHeader)
				write(DataFile, SimBetasString)
				write(DataFile, QCorrelatorsHeader)
				close(DataFile)
				
				open(FilePathOut, "a") do io
					M = QCorrelators[Scheme * "-Seq=$Sequential"]
				    writedlm(io, M, "; ")
				end
				
				# Write binned data on file
				FilePathOut = DirPathOut * "/Q_deep/$(Scheme)_NSweeps=$(NSweepsString)_QHistograms.txt"
				GeneralHeader = "# " * Scheme * ", Sequential=$(Sequential), NSweeps=" * NSweepsString * "\n"
				
				DataFile = open(FilePathOut, "w")
				write(DataFile, GeneralHeader)
				write(DataFile, QHistogramsHeader)
				close(DataFile)
				
				open(FilePathOut, "a") do io
					M = QHistograms[Scheme * "-Seq=$Sequential"]
				    writedlm(io, M, "; ")
				end
			end
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
