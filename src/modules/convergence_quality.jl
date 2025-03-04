#!/usr/bin/julia

# ------------------------------------- Run ------------------------------------

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
	MetropolisElapsedTime = 0
	
	HeatbathQMatrix = zeros(Int64, floor(Int64, NSweeps*N/QStep), length(SimBetas))
	HeatbathElapsedTime = 0
	
	for Heatbath in [false, true]
	
		ElapsedTime = @elapsed begin	# Start time tracking
		
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
		end
		
		if !Heatbath
			MetropolisQMatrix = QMatrix
			MetropolisElapsedTime = ElapsedTime
		elseif Heatbath
			HeatbathQMatrix = QMatrix
			HeatbathElapsedTime = ElapsedTime
		end
		
	end
	
	return MetropolisQMatrix, MetropolisElapsedTime, HeatbathQMatrix, HeatbathElapsedTime
end

# ----------------------------------- Modules ----------------------------------

function GetQCorrelator(k::Int64,
						QSamples::Vector{Int64})
	N = length(QSamples)
	QCorrelator = 0
	for j in 1:N-k
		QCorrelator += QSamples[j]*QSamples[j+k]
	end
	QCorrelator /= N-k
	return QCorrelator
end

function GetQBlockLengths(QSamples::Vector{Int64})
	N = length(QSamples)
	BlockLengths = []
	iStart = 1
	while iStart<N
		i = iStart+1
		while QSamples[i] == QSamples[iStart] && i<N
			i += 1
		end
		BlockLength = i-iStart	
		iStart = i+1
		append!(BlockLengths,BlockLength)
	end
	return BlockLengths
end
