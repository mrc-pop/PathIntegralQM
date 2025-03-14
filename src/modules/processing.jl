#!/usr/bin/julia

using Statistics

"""
Reshape Data by dividing them in blocks of length BlockLength, then compute the mean
value of each block and create a new matrix storing these mean values. The function
accepts both float and int values.
"""
function BlockData(Data::AbstractArray{<:Number}, BlockLength::Int64)
    # Number of blocks (using integer division ÷)
    BlockNumber = size(Data, 1) ÷ BlockLength

    # Initialize array of blocked data
    BlockedData = zeros(Float64, BlockNumber, size(Data, 2))

    for i in 1:BlockNumber
        iStart = (i - 1) * BlockLength + 1
        iEnd = i * BlockLength
        BlockedData[i, :] = mean(Data[iStart:iEnd,:], dims = 1)
    end

    return BlockedData
end

"""
Specific instance of BlockData to block Q and output Q, SigmaQ, Q², SigmaQ²
"""
function BlockQ(QQ::Array{Int64}, k::Int64)

    Q = mean(QQ)
    Q2 = mean(QQ.^2)

    BlockedQQ = BlockData(QQ, k)
    SigmaQ = std(BlockedQQ, corrected=true) / sqrt(length(BlockedQQ))

    BlockedQQ2 = BlockData(QQ.^2, k)
    SigmaQ2 = std(BlockedQQ2, corrected=true) / sqrt(length(BlockedQQ2))

    return Q, SigmaQ, Q2, SigmaQ2
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
