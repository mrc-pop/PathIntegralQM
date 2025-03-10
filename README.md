# PathIntegralQM

## Introduction

This repository contains the files needed to simulate the problem of a quantum mechanical particle on a circumference, 
using Path Integral Monte Carlo and the  Julia programming language. The work was related to the  university course
"Numerical Methods for Physics", held in AY 2023/2024 by Professor Claudio Bonati at the University of Pisa.

This work was done by Alessandro Gori ([@nepero27178](https://github.com/nepero27178)) and me.

We sample observables at finite temperature by discretizing the interval $[0,\beta]$ in $N$ imaginary time slices, and
computing the averages over all paths through a Markov-Chain Monte Carlo sampling. Every path is stored as a 
1-dimensional array. The code supports local updates (Metropolis, heatbath) as well as the "tailor" update described in 
the paper [arXiv:1709.10034](https://arxiv.org/abs/1709.10034). There is also the possibility for parallel tempering, 
which involves simulating a number $N_R$ of systems in parallel, with differrent lattice spacings, and sample a joint
probability distribution in order to help the system of interest decorrelate more quickly. The main algorithmic
subtlety is a phenomenon of a *topological critical slowing down*, due to the nontrivial topology of the configuration 
space, as well as the high energy cost of changing topological sector (i.e., change the winding number $Q$ of a path) using
local updates only.

## Outline of the repository

- The files under `/test` are used to check and debug the other codes, as well as comparing to reference benchmarks.

- The files under `/modules` contain the routines (Monte Carlo updates, plotting functions, etc.) which are used throughout
the other codes.

- The file `convergence.jl` contains an in-depth analysis of the different local updates, namely Metropolis and heatbath.
It allows to perform simulations for different lattice sizes and temperatures and analyze them through histograms of $Q$
and study the autocorrelations of the MC data quantitatively.

- The files `simulations.jl` and `simulations_parallel_tempering.jl` contain routines for simulating the system with 
different algorithms and settings. The simulation parameters are set through the `setup/simulations_setup.jl` script.

- The file `analysis.jl` is used to compute sample averages and standard deviations of quantities of interest such as
$\langle Q \rangle$ and $\langle Q^2 \rangle$, as well as plotting the results. 

