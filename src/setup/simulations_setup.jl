#!/usr/bin/julia

# TODO set simulations routine according to the best tradeoff,
# determined using the other scripts

"""
Choose model parameters; Metropolis parameters; tailor parameters for the simulation routine.
"""

const NN = [20, 40, 60, 80, 100]			# Number of lattice points
const SimBetas = [0.01, 0.025. 0.05, 0.075, 
				  0.1, 0.25, 0.5, 0.75,
				  1.0, 2.5, 5.0, 7.5,
				  10.0]          			# Adimensional inverse temperature

const Heatbath = false          			# Type of local algorithm (false → Metropolis)
const NSweepsTherm = Int(1e6)   			# Number of updates of the whole lattice for thermalization
const NSweeps = Int(1e6)        			# Number of updates of the whole lattice
const Δ = 0.5                   			# Metropolis interval width
const Sequential = true        				# Sequential or random site choice

const QSteps = fill(1,length(NN))    		# How often to compute Q, for each N (0 = never, 1=always, n=after n-1 steps)
const TailorSteps = round.(1.0 .* NN)    	# How often to propose a tailor update, for each N (0 = never)
    # fill(0,length(NN))
const ε_over_η = 0.2            			# Tolerance of tailor update in units of η
