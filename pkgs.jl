using Pkg
using CSV
using LinearAlgebra, Distributions
using JuMP
using Mosek, MosekTools
using PyCall
using Colors
using JLD
import Base: sqrt, string, abs

function string(array::Array{T, 1} where T <: Number)
    return [string(element) for element in array]
end

function sqrt(array::Vector{T} where T <: Real)
    return [sqrt(a) for a in array]
end

function abs(array::Vector{T} where T <: Number)
    return [abs(a) for a in array]
end

cd(@__DIR__)
