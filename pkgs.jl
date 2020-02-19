using Pkg
using CSV
using LinearAlgebra, Distributions
using JuMP
using Mosek, MosekTools
using PyCall
using PyPlot, Colors
using JLD
import Base: sqrt

cd(@__DIR__)
