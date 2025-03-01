#!/usr/bin/julia

# TODO set simulations routine according to the best tradeoff,
# determined using the other scripts

"""
Choose model parameters; metropolis parameters; tailor parameters for the simulation routine.
"""

const NN = [20]                   # number of lattice points
const SimBetas = [1.0, 2.0]       # adimensional inverse temperature

const NSweepsTherm = Int(1e4)   # number of updates of the whole lattice for thermalization
const NSweeps = Int(1e8)        # number of updates of the whole lattice
const Δ = 0.5                   # metro interval width
const sequential = false        # sequential or random site choice

const TailorStep = 50           # how often to perform tailor update (0 = never)
const ε_over_η = 1.0            # tolerance of tailor update in units of η
