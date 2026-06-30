using MutationLoadDynamics
using Test
using Random
using StatsBase
using AbstractTrees
using Distributions
using Statistics
using DataStructures

tests = [
    "initialisation",
    "events",
    "simulations",
    "statistics",
]

@testset "MutationLoadDynamics.jl" begin
    for test in tests
        @testset "$test" begin
            include(test * ".jl")
        end
    end
end
