"""
    initialize_population(; fitness_init=1.0, time=0.0) -> Population

Create a population containing a single founding cell with the given initial fitness.
The cell acquires no driver mutations at birth (it is the root of the lineage tree).

```julia
pop = initialize_population(fitness_init = 1.0)
```
"""
function initialize_population(; fitness_init::Float64 = 1.0, time::Float64 = 0.0)
    cell = NonMarkovCell(1, time, 0, fitness_init)
    node = BinaryNode(cell)
    cells = Dict{Int64, BinaryNode{NonMarkovCell}}(1 => node)
    return Population(cells, time, 1)
end

"""
    initialize_population(N::Int; fitness_init=1.0, time=0.0) -> Population

Create a population of `N` identical independent cells, each with `fitness_init`.
Useful for starting from a pre-existing pool rather than a single founder.
"""
function initialize_population(N::Int; fitness_init::Float64 = 1.0, time::Float64 = 0.0)
    cells = Dict{Int64, BinaryNode{NonMarkovCell}}()
    for id in 1:N
        cell = NonMarkovCell(id, time, 0, fitness_init)
        cells[id] = BinaryNode(cell)
    end
    return Population(cells, time, N)
end
