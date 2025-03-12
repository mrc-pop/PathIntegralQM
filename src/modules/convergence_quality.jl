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

	for Heatbath in [false]#, true]

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

# ----------------------------- In-depth analysis ------------------------------

function RunDeepAnalysis(
	DirPathOut::String,
	NSweepsString::String,
	InputMetroData::Matrix{Int64},
	InputHeatData::Matrix{Int64},
	SimBetas::Vector{Float64},
	Sequential::Bool,
	AdditionalHeaders::Dict{String, String};
	kMax = 100,
	Skip = 5,
	Bins = -0.5:1.0:100.5 # 1, ..., 99, 100
)

	RowsNumber = floor(Int64, (kMax+1)/Skip)
	MetropolisQCorrelators = zeros(RowsNumber,length(SimBetas))
	HeatbathQCorrelators = zeros(RowsNumber,length(SimBetas))
	QCorrelators = Dict()

	MetropolisQBlockLengths = Any[]
	HeatbathQBlockLengths = Any[]

	MetropolisHistogram = Any[]
	HeatbathHistogram = Any[]
	QHistograms = Dict()

	for (sb,SimBeta) in enumerate(SimBetas)

		TmpCounter = 1
		for k in 0:Skip:kMax
			MetropolisQCorrelators[TmpCounter,sb] = GetQCorrelator(k, InputMetroData[:,sb])
			HeatbathQCorrelators[TmpCounter,sb] = GetQCorrelator(k, InputHeatData[:,sb])
			TmpCounter += 1
		end

		MetroQBL = GetQBlockLengths(InputMetroData[:,sb])
		push!(MetropolisQBlockLengths, MetroQBL)

		HeatQBL = GetQBlockLengths(InputHeatData[:,sb])
		push!(HeatbathQBlockLengths, HeatQBL)

		hM = fit(Histogram, MetroQBL, Bins)
		push!(MetropolisHistogram, [SimBeta, hM.weights, hM.edges[1]])

		hH = fit(Histogram, HeatQBL, Bins)
		push!(HeatbathHistogram, [SimBeta, hH.weights, hH.edges[1]])

	end

	QCorrelators = Dict([("Metropolis-Seq=$Sequential", MetropolisQCorrelators),
						("Heatbath-Seq=$Sequential", HeatbathQCorrelators)])

	QHistograms = Dict([("Metropolis-Seq=$Sequential", MetropolisHistogram),
						 ("Heatbath-Seq=$Sequential", HeatbathHistogram)])

	for (s,Scheme) in enumerate(["Metropolis", "Heatbath"])
		mkpath(DirPathOut * "/Q_deep/")

		# Write correlators on file
		FilePathOut = DirPathOut * "/Q_deep/$(Scheme)_NSweeps=$(NSweepsString)_QCorrelators.txt"
		GeneralHeader = "# " * Scheme * ", Sequential=$(Sequential), NSweeps=" * NSweepsString * ", Skip=$(Skip)\n"

		DataFile = open(FilePathOut, "w")
		write(DataFile, GeneralHeader)
		write(DataFile, AdditionalHeaders["SimBetas"])
		write(DataFile, AdditionalHeaders["SimBetasValues"])
		write(DataFile, AdditionalHeaders["QCorrelators"])
		close(DataFile)

		open(FilePathOut, "a") do io
			M = QCorrelators[Scheme * "-Seq=$Sequential"]
		    writedlm(io, M, "; ")
		end

		# Write binned data on file
		FilePathOut = DirPathOut * "/Q_deep/$(Scheme)_NSweeps=" * NSweepsString * "_QHistograms.txt"
		GeneralHeader = "# " * Scheme * ", Sequential=$(Sequential), NSweeps=" * NSweepsString * "\n"

		DataFile = open(FilePathOut, "w")
		write(DataFile, GeneralHeader)
		write(DataFile, AdditionalHeaders["QHistograms"])
		close(DataFile)

		open(FilePathOut, "a") do io
			M = QHistograms[Scheme * "-Seq=$Sequential"]
		    writedlm(io, M, "; ")
		end
	end
end

# ----------------------------------- Modules ----------------------------------

function GetQCorrelator(k::Int64,
						QSamples::Vector{Int64})
	N = length(QSamples)
	AvgQ = mean(QSamples)
	StdQ = std(QSamples)
	QCorrelator = 0
	for j in 1:N-k
		QCorrelator += (QSamples[j]-AvgQ) * (QSamples[j+k]-AvgQ)
	end
	
	if StdQ != 0
		QCorrelator /= (StdQ^2 * (N-k))
	end
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

# ---------------------------- Q  Path processing ------------------------------

function ProcessQDataFilePath(
	FilePathIn::String
)

	# Find directory delimiters "/"
	DirDelims = findall('/', FilePathIn)

	# Check if the files was generated sequentially or randomly
	Sequential = missing
	SequentialString = FilePathIn[ DirDelims[end-2]+1 : DirDelims[end-1]-1 ]
	if SequentialString == "sequential"
		Sequential = true
	elseif SequentialString == "random"
		Sequential = false
	else
		error("Invalid FilePathIn. Impossible to distinguish sequential/random.")
	end

	# Extract N reading the directory name
	NString = FilePathIn[ DirDelims[end-1]+1 : DirDelims[end]-1 ]
	N = parse(Int64, NString[3:end])

	# Extract the scheme reading the file name
	FileString = FilePathIn[ DirDelims[end]+1 : end ]
	UnderScoreDelims = findall('_', FileString)
	Scheme = FileString[ 1 : UnderScoreDelims[1]-1 ]

	# Extract NSweeps reading the file name
	if length(UnderScoreDelims)==1
		NSweepsString = FileString[ UnderScoreDelims[1]+9 : end-4 ]
		NSweeps = Int64(parse(Float64, NSweepsString))
	elseif length(UnderScoreDelims)==2
		NSweepsString = FileString[ UnderScoreDelims[1]+9 : UnderScoreDelims[2]-1 ]
		NSweeps = Int64(parse(Float64, NSweepsString))
	end

	return Sequential, N, Scheme, NSweeps

end
