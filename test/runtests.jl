using MutationLoadDynamics
using Test
using Random
using AbstractTrees
using Distributions
using Statistics
using DataStructures
using StableRNGs

include("fixtures.jl")

tests = [
    "initialisation",
    "events",
    "simulations",
    "regression",
    "chaining",
    "measurements",
    "statistics",
    "validation",
    "sampling",
]

@testset "MutationLoadDynamics.jl" begin
    for test in tests
        @testset "$test" begin
            include(test * ".jl")
        end
    end
end
