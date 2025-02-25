#!/usr/bin/julia

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

function main()
	
	NSteps = 100000						# Number of single-site update steps
    N = 20								# Number of time steps
	SimBeta = 1.0						# Adimensional beta (Beta/kB*T)

    Config = SetLattice(SimBeta, N)		# Initalize path
	IgnoredSteps = 10
	WindingNumbers = zeros(Int64, floor(Int64, NSteps/IgnoredSteps), 2)
	IgnoreCounter = IgnoredSteps

	Schemes = ["Metropolis", "Heatbath"]
	
	Metropolis = false 					# TODO Change
	AccSteps = 0    

	for (j,Scheme) in enumerate(Schemes)
		println("Performing $NSteps $Scheme steps...")

		for i in 1:NSteps
		    Site = i % (N-2) + 2 # from 2 to N-1 (sequential)
		    # Site = rand(2:N-1) # (random)
		    
		    if Scheme=="Metropolis"
		    	AccSteps += MetropolisUpdate!(Config, Site; Δ=0.05, verbose=false)
	 		elseif Scheme=="Heatbath"
	 			HeatBathUpdate!(Config, Site; verbose=false)
	 		else
	 			println("Invalid algorithm. Exiting.")
	 			exit()
	 		end		
	 		       
		    if IgnoreCounter==0 
		   		IgnoreCounter = IgnoredSteps
		   		WindingNumbers[ceil(Int64, i/IgnoredSteps), j] = round(Int64, CalculateQ(Config))
		   	elseif IgnoreCounter>0
		   		IgnoreCounter -= 1
		   	else
		   		println("There's something strange with your counter, man.")
		   	end
		end

		if Metropolis
			println("Accepted steps: $AccSteps/$NSteps")
		end

	end
	
	Waiting=true
	print("Plot results? (y/n) ")
	UserPlot = readline()
	
	while Waiting
		if UserPlot=="y"
			Waiting=false
			# TODO Sistema
			QHistogram = histogram(
				[WindingNumbers[:,1] WindingNumbers[:,2]],
				bins=range(-4.5,4.5,length=10),
				ylims=[0,1],
				normalize=:pdf,
				color=[MyColors[1] MyColors[2]],
				xlabel=L"$Q$ (winding number)",
				ylabel="Normalized occurrencies",
				title=L"Single-site full-chain sequential updates ($\tilde{\beta}=%$(SimBeta)$, $N=%$N$)",
				label=["Metropolis" "Heatbath"]
			)
			savefig(QHistogram, PROJECT_ROOT * "/../convergence/QHistogram_SimBeta=$(SimBeta)_N=$N.pdf")
			
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
