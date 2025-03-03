#!/usr/bin/julia

using DelimitedFiles
using Dates

PROJECT_ROOT = @__DIR__
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
	Scheme::String,
	UserDelta::Float64,
	N::Int64,
	NSteps::Int64,
	SimBetas::Vector{Float64},
	IgnoredSteps::Int64;
	RandomSite=false
)
	
	QScheme = zeros(Int64, floor(Int64, NSteps/(IgnoredSteps+1)), length(SimBetas))
	
	# Run simulations
	for	(i,SimBeta) in enumerate(SimBetas)
	
		QCounter = 1
		
		Config = SetLattice(SimBeta, N)	# Initalize path
		printstyled("Performing simulation $i/$(length(SimBetas)) at SimBeta=$(SimBetas[i])...", color=:yellow)

		IgnoreCounter = IgnoredSteps+1	# Initialize ignore counter
		for j in 1:NSteps
		
			if !RandomSite
				Site = j % (N-3) + 2 	# from 2 to N-1 (sequential)
			elseif RandomSite
				Site = rand(2:N-1) 		# (random)
			end
			
			# TODO Vary UserDelta
			
			if Scheme=="Metropolis"
				MetropolisUpdate!(Config, Site; Δ=UserDelta, verbose=false)
			elseif Scheme=="Heatbath"
				HeatBathUpdate!(Config, Site; verbose=false)
			else
				error("Invalid algorithm. Exiting.")
			end
	 		
			if IgnoreCounter==0 
		   		IgnoreCounter = IgnoredSteps
		   		QScheme[QCounter, i] = round(Int64, CalculateQ(Config))
		   		QCounter += 1
		   		
		   	elseif IgnoreCounter>0
		   		IgnoreCounter -= 1
		   	else
		   		printstyled("There's something strange with your counter, man.\n", color=:red)
		   	end	
		end
		
		printstyled(" Done.\n", color=:green)
	end

	return QScheme
end

function PlotHistograms(
	DirPathIn::String,
	Scheme::String,
	UserDelta::Float64,
	N::Int64,
	NSteps::Int64,
	SimBetas::Vector{Float64},
	IgnoredSteps::Int64
)

	FilePathIn = DirPathIn * "/$(Scheme)_NSteps=$(NSteps).txt"
	QData = readdlm(FilePathIn, ';', comments=true)
	
	# Mastruzzo here
	PlotDict = Vector{Dict}()
	
	for (j,SimBeta) in enumerate(SimBetas)
	
		h = histogram(
			QData[:,j],
			normalize=:pdf,
			fillcolor=MyColors[1],
			label=L"$\tilde{\beta}=%$SimBeta$",
			#
			bins=range(-4.5,4.5,length=10),
			ylims=[0,1],
			size=(440,300),
			title=L"%$Scheme, sequential ($\eta=%$(SimBetas[j]/N)$)",
			xlabel=L"$Q$ (winding number)",
			ylabel="Normalized occurrencies",
		)
		
		FilePathOut = DirPathIn * "$(Scheme)_SimBeta=$SimBeta.pdf"
		savefig(h, FilePathOut)
		
		push!(PlotDict, Dict("Histogram"=>h))
		
	end
end

function main()

	RS = false												# Random site selection
	
	NStepsString = "1E5"									# Number of single-site update steps
	NSteps = Int64(parse(Float64, NStepsString))
	N = 20													# Number of time steps
	SimBetas = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]		# Adimensional betas (Beta/kB*T)
	# Etas = round.(SimBetas./N, digits=2)

	# Write files headers
	SimBetasHeader = "# First line: SimBetas\n"
	QHeader = "# Following lines: Qs (divided for each eta) [calculated $(now())]\n"
	
	SimBetasString = ""
	for SimBeta in SimBetas
		SimBetasString *= "$SimBeta; "
	end
	SimBetasChars = collect(SimBetasString)
	popat!(SimBetasChars, length(SimBetasString)-1)
	SimBetasString = String(SimBetasChars) * "\n"

	IgnoredSteps = 0						# After 50 steps extract Q

	Schemes = ["Metropolis", "Heatbath"] 	# Algorithms
	UserDelta = 0.05

	if RS
		DirPathOut = PROJECT_ROOT * "/../convergence/data/random/N=$(N)/"
	elseif !RS
		DirPathOut = PROJECT_ROOT * "/../convergence/data/sequential/N=$(N)/"
	end
	mkpath(DirPathOut)

	for (s,Scheme) in enumerate(Schemes)
		println("Starting $(Scheme) simulations...")

		FilePathOut = DirPathOut * "/$(Scheme)_NSteps=$(NStepsString).txt"
		DataFile = open(FilePathOut, "w")
		write(DataFile, SimBetasHeader)
		write(DataFile, SimBetasString)
		write(DataFile, QHeader)
		close(DataFile)
		
		QScheme = RunConvergenceSimulations(Scheme, UserDelta, N, NSteps, SimBetas, IgnoredSteps; RandomSite=RS)
		
		# Write on file
		DataFile = open(FilePathOut, "a")
		
		for i in 1:floor(Int64, NSteps/(IgnoredSteps+1))
			# Generate entry
			Entry = ""
		
			for j in 1:length(SimBetas)
				Entry *= "$(QScheme[i,j]); "
			end
		
			EntryChars = collect(Entry)
			pop!(EntryChars)
			EntryChars[end] = '\n'
			Entry = String(EntryChars)
			write(DataFile, Entry)
		end
		
		close(DataFile)
	end
	
	Waiting=true
	print("Plot results? (y/n) ")
	UserPlot = readline()
	
	while Waiting
		if UserPlot=="y"
			Waiting=false
			for Scheme in Schemes
				PlotHistograms(DirPathOut, Scheme, UserDelta, N, NSteps, SimBetas, IgnoredSteps)
			end
			
		elseif UserPlot=="n"
			Waiting=false
			println("Aborting.")
			exit()
		else
			Waiting=true
			print("Invalid input. Plotting results? (y/n) ")
			UserPlot = readline()
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
