#!/usr/bin/julia

using DelimitedFiles
using Statistics

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../src/modules/processing.jl")

const FilePath = PROJECT_ROOT*"/../simulations/data_Q_mock.txt"
const k = 10000  # block length
const R = 100    # resamples

function main()

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
