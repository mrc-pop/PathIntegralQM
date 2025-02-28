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
	
		Config = SetLattice(SimBeta, N)	# Initalize path
		println("Performing simulation $i/$(length(SimBetas)): $NSteps $Scheme steps at Eta=$(SimBetas[i]/N)...")
		
		IgnoreCounter = IgnoredSteps+1	# Initialize ignore counter
		for j in 0:(NSteps-1)
		
			if !RandomSite
				Site = j % (N-3) + 2 	# from 2 to N-1 (sequential)
			elseif RandomSite
				Site = rand(2:N-1) 	# (random)
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
		   		IgnoreCounter = IgnoredSteps+1
		   		QScheme[ceil(Int64, j/(IgnoredSteps+1)), i] = round(Int64, CalculateQ(Config))
		   
		   	elseif IgnoreCounter>0
		   		IgnoreCounter -= 1
		   	else
		   		println("There's something strange with your counter, man.")
		   	end	
		end
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
	
	# display(PlotDict)
	
#	MetropolisHistogram = plot(
#		PlotDict[1]["Histogram"],
#		PlotDict[2]["Histogram"],
#		PlotDict[3]["Histogram"],
#		PlotDict[4]["Histogram"],
#		PlotDict[5]["Histogram"],
#		layout=(5,1),
##		xlabel=L"$Q$ (winding number)",
##		ylabel="Normalized occurrencies",
#	)
#	FilePathOut = DirPathIn * "$(Scheme).pdf"
#	savefig(MetropolisHistogram, FilePathOut)
end

function main()

	RS = true								# Random site selection
	
	NSteps = Int64(1E8)						# Number of single-site update steps
    N = 100									# Number of time steps
	SimBetas = [0.1, 0.5, 1.0, 2.0, 4.0]	# Adimensional betas (Beta/kB*T)
	Etas = SimBetas./N

	# Write files headers
	Header = "# "
	for Eta in Etas
		Header *= "Q-Eta=$Eta; "
	end
	HeaderChars = collect(Header)
	popat!(HeaderChars, length(HeaderChars)-1)
	Header = String(HeaderChars)

	IgnoredSteps = 49						# After 50 steps extract Q

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

		FilePathOut = DirPathOut * "/$(Scheme)_NSteps=$(NSteps).txt"
		DataFile = open(FilePathOut, "w")
		write(DataFile, Header * "[calculated $(now())]\n")
		close(DataFile)
		
		QScheme = RunConvergenceSimulations(Scheme, UserDelta, N, NSteps, SimBetas, IgnoredSteps; RandomSite=RS)
		
		# Write on file
		DataFile = open(FilePathOut, "a")
		
		for i in 1:floor(Int64, NSteps/(IgnoredSteps+1))
			# Generate entry
			Entry = ""
		
			for j in 1:length(Etas)
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
