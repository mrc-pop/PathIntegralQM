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

function main()

	BiasSimBeta = 2.0
	BiasN = 50
	SequentialRatios = [1, 2, 5, 10]
	@info "Model settings" BiasN BiasSimBeta SequentialRatios

	NSweepsString = @sprintf "%.1e" NSweeps
	DirPathOut = PROJECT_ROOT * "/../convergence/sequential/biased"
	mkpath(DirPathOut)
	
	BiasDirPathOut = DirPathOut

	for (n,N) in enumerate([BiasN])

		for (sr,SequentialRatio) in enumerate(SequentialRatios)
		
			QStep = 1
			DirPathOut = BiasDirPathOut * "/N=$(BiasN)/"
			mkpath(DirPathOut)

			println()
			@info "Starting biased simulations for N=$(BiasN)..."
			Results = RunBiasedSequentialSimulations(BiasN, NSweepsTherm, NSweeps, [BiasSimBeta], QStep, SequentialRatio)
			MetropolisQMatrix, HeatbathQMatrix = Results
			QMatrices = Dict([("Metropolis", MetropolisQMatrix),
							  ("Heatbath", HeatbathQMatrix)])

			for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])

				@info "Writing simulation $(2*(sr-1)+s)/$(length(SequentialRatios)*2)" Scheme SequentialRatio

				FilePathOut = DirPathOut *
					"/" * Scheme * 
					"_SequentialRatio=$(SequentialRatio)" *
					"_NSweeps=" * NSweepsString * ".txt"
					
				GeneralHeader = "# " * Scheme *
								", NSweeps=" * NSweepsString *
								", SequentialRatio=$(SequentialRatio)\n"

				DataFile = open(FilePathOut, "w")
				write(DataFile, GeneralHeader)
				close(DataFile)

				open(FilePathOut, "a") do io
				    writedlm(io, QMatrices[Scheme], "; ")
				end
			end

			# -------------------------- In-depth analysis ---------------------

			kMax = 100
			Skip = 5

			RowsNumber = floor(Int64, kMax/Skip)+1
			MetropolisQCorrelators = zeros(RowsNumber,length([BiasSimBeta]))
			HeatbathQCorrelators = zeros(RowsNumber,length([BiasSimBeta]))
			QCorrelators = Dict()

			for (sb,SimBeta) in enumerate([BiasSimBeta])

				TmpCounter = 1
				for k in 0:Skip:kMax
					MetropolisQCorrelators[TmpCounter,sb] = GetQCorrelator(k, MetropolisQMatrix[:,sb])
					HeatbathQCorrelators[TmpCounter,sb] = GetQCorrelator(k, HeatbathQMatrix[:,sb])
					TmpCounter += 1
				end
				
			end
			
			QCorrelators = Dict([("Metropolis", MetropolisQCorrelators),
								("Heatbath", HeatbathQCorrelators)])
			
			for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])

				FilePathOut = DirPathOut *
					"/" * Scheme *
					"_SequentialRatio=$(SequentialRatio)" *
					"_NSweeps=" * NSweepsString *
					"_QCorrelators.txt"
				GeneralHeader = "# " * Scheme *
								", Sequential=true," *
								" NSweeps=" * NSweepsString * "\n"

				DataFile = open(FilePathOut, "w")
				write(DataFile, GeneralHeader)
				write(DataFile, "# SimBeta\n")
				write(DataFile, "$(BiasSimBeta)\n")
				write(DataFile, "# Q Correlators\n")
				close(DataFile)

				open(FilePathOut, "a") do io
					M = QCorrelators[Scheme]
				    writedlm(io, M, "; ")
				end
			end
		end
	end

	# ------------------------------ Detailed plot -----------------------------

	pgfplotsx()
	PlotBiasedQCorrelatorsCompared(DirPathOut, BiasSimBeta, BiasN, NSweeps, SequentialRatios)

end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
