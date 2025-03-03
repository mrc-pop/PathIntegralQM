#!/usr/bin/julia

using Printf
using DelimitedFiles
using Dates

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/simulations_setup.jl")
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/modules/montecarlo.jl")
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

function RunConvergenceSimulations(
	N::Int64,
	NSweepsTherm::Int64,
	NSweeps::Int64,
	SimBetas::Vector{Float64},
	QStep::Int64;
	Sequential=false,
	UserDelta=Δ
)

	MetropolisQMatrix = zeros(Int64, floor(Int64, NSweeps*N/QStep), length(SimBetas))
	HeatbathQMatrix = MetropolisQMatrix
	
	for Heatbath in [false, true]
		
		if !Heatbath
			@info "Starting Metropolis simulations..."
		elseif Heatbath
			@info "Starting Heatbath simulations..."
		end
	
		QMatrix = zeros(Int64, floor(Int64, NSweeps*N/(QStep)), length(SimBetas))
	
		if !Heatbath
			QMatrix = MetropolisQMatrix
		elseif Heatbath
			QMatrix = HeatbathQMatrix
		end
		
		Scheme = Heatbath ? "Heatbath" : "Metropolis"
	
		# Run simulations
		for	(sb,SimBeta) in enumerate(SimBetas)
			
			QCounter = 1
			Config = SetLattice(SimBeta, N)	# Initalize path
			
			println()
			@info "Model settings" N SimBeta
		    @info Scheme * " settings" Heatbath Δ Sequential NSweepsTherm NSweeps
			@info "Simulation $sb/$(length(SimBetas))"

			# Thermalization
		    println("\nPerforming $NSweepsTherm " * Scheme * " sweeps for thermalization...")
		    @time for i in 1:NSweepsTherm
		        CurrentSweepSteps = N*(i-1)
		        for j in 1:N
		            # Choose site
		            if Sequential
		                Site = mod1(CurrentSweepSteps + j, N)
		            elseif !Sequential
		                Site = rand(1:N)
		            end
		            # Perform update
		            if !Heatbath
		            	MetropolisUpdate!(Config, Site; Δ=UserDelta)
		            elseif Heatbath
		            	HeatBathUpdate!(Config, Site)
		            end         
		        end
		    end
		    
	#        unicodeplots()
	#        p = PlotPathUnicode(Config)
	#        @info "Plot after thermalization " p

		    # Local update sweeps
		    println("\nPerforming $NSweeps "* Scheme * " sweeps of the whole lattice...")

			# Main run
			@time for i in 1:NSweeps
				CurrentSweepSteps = N*(i-1)
				for j in 1:N
					# Choose site
				    if Sequential
				        Site = mod1(CurrentSweepSteps + j, N)
				    elseif !Sequential
				        Site = rand(1:N)
				    end
					# Perform update
		            if !Heatbath
		            	MetropolisUpdate!(Config, Site; Δ=UserDelta) 
		            elseif Heatbath
		            	HeatBathUpdate!(Config, Site)
		            end
			 		# Measure Q
				    if QStep !== 0 && mod(CurrentSweepSteps + j, QStep) == 0
				        QMatrix[QCounter,sb] = CalculateQ(Config)
				        QCounter += 1
				    end
				end
			end
		end
		
		if !Heatbath
			MetropolisQMatrix = QMatrix
		elseif Heatbath
			HeatbathQMatrix = QMatrix
		end
		
	end
	
	return MetropolisQMatrix, HeatbathQMatrix
end

# ------------------------------------ Main ------------------------------------ 

function main()

	# Prepare files headers: first line contains SimBetas
	SimBetasHeader = "# First line: SimBetas\n"
	QHeader = "# Following lines: Qs (divided for each eta) [calculated $(now())]\n"
	
	SimBetasString = ""
	for SimBeta in SimBetas
		SimBetasString *= "$SimBeta; "
	end
	SimBetasChars = collect(SimBetasString)
	popat!(SimBetasChars, length(SimBetasString)-1)
	SimBetasString = String(SimBetasChars) * "\n"

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
		MetropolisQMatrix, HeatbathQMatrix = RunConvergenceSimulations(N, NSweepsTherm, NSweeps, SimBetas, QStep; Sequential)
		QMatrices = Dict([("Metropolis-Seq=$Sequential", MetropolisQMatrix), 
						  ("Heatbath-Seq=$Sequential", HeatbathQMatrix)])

		for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])

			NSweepsString = @sprintf "%.1e" NSweeps
			FilePathOut = DirPathOut * "/$(Scheme)_NSweeps=$(NSweepsString).txt"
			GeneralHeader = "# " * Scheme * ", Sequential=$(Sequential), NSweeps=" * NSweepsString * "\n"
			
			DataFile = open(FilePathOut, "w")
			write(DataFile, GeneralHeader)
			write(DataFile, SimBetasHeader)
			write(DataFile, SimBetasString)
			write(DataFile, QHeader)
			close(DataFile)
			
			open(FilePathOut, "a") do io
				A = QMatrices[Scheme * "-Seq=$Sequential"]
		        writedlm(io, A, "; ")
		    end
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
