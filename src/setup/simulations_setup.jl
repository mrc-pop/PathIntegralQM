#!/usr/bin/julia

# TODO set simulations routine according to the best tradeoff,
# determined using the other scripts

"""
Choose model parameters; metropolis parameters; tailor parameters for the simulation routine.
"""

const NN = [20, 40, 60, 80, 100, 200, 300]              # number of lattice points
const SimBetas = [2.0]          # adimensional inverse temperature

const heatbath = false          # type of local algorithm (false → metropolis)
const NSweepsTherm = Int(1e6)   # number of updates of the whole lattice for thermalization
const NSweeps = Int(1e6)        # number of updates of the whole lattice
const Δ = 0.5                   # metro interval width
const sequential = true        # sequential or random site choice

const TailorSteps = round.(1.0 .* NN)    # how often to propose tailor update, for each N (0 = never)
    # fill(0,length(NN))
const ε_over_η = 0.2            # tolerance of tailor update in units of η
