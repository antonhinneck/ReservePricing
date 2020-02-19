using Pkg
using CSV
using LinearAlgebra, Distributions
using JuMP
using Mosek, MosekTools
using PyCall
using PyPlot, Colors
using JLD
import Base: sqrt, string

function string(array::Array{T, 1} where T <: Number)
    return [string(element) for element in array]
end

cd(@__DIR__)
