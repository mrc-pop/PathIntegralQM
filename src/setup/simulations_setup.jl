#!/usr/bin/julia

# TODO set simulations routine according to the best tradeoff,
# determined using the other scripts

"""MODEL SETTINGS"""
const NN = [50, 400]
const SimBetas = vcat(0.1, [x for x in 0.5:0.5:10.0])

"""
const NN = [25, 50, 75, 100,
            125, 150, 175, 200,
            300, 400]			            # Number of lattice points

const SimBetas = [0.01, 0.025, 0.05, 0.075,
				  0.1, 0.25, 0.5, 0.75,
				  1.0, 2.5, 5.0, 7.5,
				  10.0]          			# Adimensional inverse temperature
"""

"""ALGORITHM SETTINGS"""
const Heatbath = false          			# Type of local algorithm (false → Metropolis)
const NSweepsTherm = Int(1e2)   			# Number of updates of the whole lattice for thermalization
const NSweeps = Int(5e5)        			# Number of updates of the whole lattice
const Δ = 0.5                   			# Metropolis interval width
const Sequential = true        				# Sequential or random site choice

const TailorSteps = round.(0.0 .* NN)       # How often to propose a tailor update, for each N (0 = never)
const ε_over_η = 0.2            			# Tolerance of tailor update in units of η

"""MEASUREMENT SETTINGS"""
const QSteps = 10 .* NN # fill(,length(NN)) # How often to compute Q, for each N (0 = never, 1=always, n=after n-1 steps)

"""PARALLEL TEMPERING SETTINGS"""
const NR = 1                                # Number of systems to be simulated in parallel
const Ratio = 1.2                           # Ratio between the η of successive systems.
const SwapStep = 50                         # How often to propose an exchange
