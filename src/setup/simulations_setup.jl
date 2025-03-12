#!/usr/bin/julia

# TODO set simulations routine according to the best tradeoff,
# determined using the other scripts

"""CONVERGENCE SETTINGS"""
const ConvergenceNNDeep = [400] # [50, 400]
const ConvergenceNNExtended = [50, 400] # [x for x in 50:50:400]
const ConvergenceSimBetasDeep = [0.5, 10.0] # [0.5, 10.0]
const ConvergenceSimBetasExtended = vcat(0.1, [x for x in 0.5:0.5:10.0])

"""MODEL SETTINGS"""
const NN = [75, 100, 125, 150, 175, 200, 300, 400]			            # Number of lattice points

const SimBetas = [2.0]        			# Adimensional inverse temperature

"""ALGORITHM SETTINGS"""
const Heatbath = false          			# Type of local algorithm (false → Metropolis)
const NSweepsTherm = 0 # Int(1e5)   			# Number of updates of the whole lattice for thermalization
const NSweeps = Int(1e7) # Int(1e6)        			# Number of updates of the whole lattice
const Δ = 0.5                   			# Metropolis interval width
const Sequential = false        				# Sequential or random site choice

const TailorSteps = round.(5 .* NN)       # How often to propose a tailor update, for each N (0 = never)
const ε_over_η = 1.0            			# Tolerance of tailor update in units of η

"""MEASUREMENT SETTINGS"""
const QSteps = NN		# How often to compute Q, for each N (0 = never, 1=always, n=after n-1 steps)
const MeasureEvery = 10 .* NN

"""PARALLEL TEMPERING SETTINGS"""
const NR = 1                                # Number of systems to be simulated in parallel
const Ratio = 1.2                           # Ratio between the η of successive systems.
const SwapStep = 50                         # How often to propose an exchange
