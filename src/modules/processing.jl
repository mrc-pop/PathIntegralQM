#!/usr/bin/julia

using DelimitedFiles
using Statistics

"""
Reshape Data by dividing them in blocks of length BlockLength, then compute the mean
value of each block and create a new matrix storing these mean values.
"""
function BlockData(Data::Matrix{Float64}, BlockLength::Int64)
    # Number of blocks
    BlockNumber = Int(floor(size(Data, 1) / BlockLength))

    # Initialize array of blocked data
    BlockedData = zeros(BlockNumber, size(Data, 2))

    for i in 1:BlockNumber
        iStart = (i - 1) * BlockLength + 1
        iEnd = i * BlockLength
        BlockedData[i, :] = mean(Data[iStart:iEnd,:], dims = 1)
    end

    return BlockedData
end

"""
Resample data for using the bootstrap method: construct a new dataset, extracting the
available data at random, with replacement.
"""
function ResampleBootstrap(Data::Matrix{Float64})
    FakeData = zeros(size(Data))

    for i in 1:size(Data,1)
        RandomIndex = rand(1:size(Data,1))
        FakeData[i,:] = Data[RandomIndex,:]
    end

    return FakeData
end

"""
Calculate the errors on Data making R bootstrap resamples and taking their std.
"""
function GetBootstrapErrors(Data::Matrix{Float64}, R::Int64)
	# Get errors using bootstrap algorithm
	FakeObservables = zeros(R, size(Data,2))

	for i in 1:R
		# Resample data to produce matrix of fake data
		FakeData = ResampleBootstrap(Data)
        # TODO if we need to manipulate the resampled data, we can do it here
        # (GetSecondaryObservables)
	end

    # Compute std dev over direction 1 (y direction)
	BootstrapErrors = std(FakeData,dims=1)

	return BootstrapErrors
end

function Jackknife()
    # TODO
end

function main()
    # SETTINGS
    FilePath=PROJECT_ROOT*"/../../simulations/data_Q_mock.txt"
    k = 10000 # block length
    R = 100 # resamples
    # end settings

    Data = readdlm(FilePath, ',', '\n'; comments=true)

    println("Reading n=$(length(Data)) values of Q. Using blocks of length $k.")

    QQ = Data
    QQ2 = QQ.^2

    println("⟨Q⟩ = ", mean(QQ))
    println("⟨Q²⟩ = ", mean(QQ2))

    BlockedQQ = BlockData(QQ, k)
    BlockedQQ2 = BlockData(QQ2, k)

    SigmaQ = std(BlockedQQ, corrected=true) / sqrt(length(BlockedQQ))
    SigmaQ2 = std(BlockedQQ2, corrected=true) / sqrt(length(BlockedQQ2))

    println("δQ = ", SigmaQ)
    println("δQ² = ", SigmaQ2)
end

main()
