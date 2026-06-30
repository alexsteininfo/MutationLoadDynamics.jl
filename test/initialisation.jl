@testset "single-cell initialization" begin
    pop = initialize_population(fitness_init = 2.0)
    @test popsize(pop) == 1
    @test pop.t == 0.0
    @test pop._next_id == 1
    cell = first(values(pop.cells)).data
    @test cell.id == 1
    @test cell.fitness ≈ 2.0
    @test cell.mutations == 0
    @test cell.birthtime ≈ 0.0
end

@testset "N-cell initialization" begin
    pop = initialize_population(5; fitness_init = 1.5)
    @test popsize(pop) == 5
    @test pop._next_id == 5
    for node in values(pop.cells)
        @test node.data.fitness ≈ 1.5
        @test node.data.mutations == 0
    end
end

@testset "allcells length and type" begin
    pop = initialize_population(3)
    cells = allcells(pop)
    @test length(cells) == 3
    @test eltype(cells) == BinaryNode{NonMarkovCell}
end

@testset "default fitness is 1.0" begin
    pop = initialize_population()
    @test first(values(pop.cells)).data.fitness ≈ 1.0
end
