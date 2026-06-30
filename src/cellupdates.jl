"""
    celldivision!(population, parent_node, t, block, rng) -> (d1, d2)

Replace `parent_node` (which has just divided) with two daughter cells in the tree and
in `population.cells`. Each daughter independently draws `j ~ Poisson(ν)` driver
mutations; for each mutation, `δ ~ block.driver_dist` and fitness is updated via
`block.fitness_update`. Returns the two daughter `BinaryNode`s.
"""
function celldivision!(
    population::Population,
    parent_node::BinaryNode{NonMarkovCell},
    t::Float64,
    block::NonMarkovBlock,
    rng::AbstractRNG,
)
    parent_fitness = parent_node.data.fitness

    d1_data = _make_daughter(population, t, parent_fitness, block, rng)
    d2_data = _make_daughter(population, t, parent_fitness, block, rng)

    d1_node = leftchild!(parent_node, d1_data)
    d2_node = rightchild!(parent_node, d2_data)

    delete!(population.cells, parent_node.data.id)
    population.cells[d1_node.data.id] = d1_node
    population.cells[d2_node.data.id] = d2_node

    return d1_node, d2_node
end

function _make_daughter(
    pop::Population,
    t::Float64,
    parent_fitness::Float64,
    block::NonMarkovBlock,
    rng::AbstractRNG,
)
    j = rand(rng, Poisson(block.ν))
    f = parent_fitness
    for _ in 1:j
        δ = rand(rng, block.driver_dist)
        f = block.fitness_update(f, δ)
    end
    pop._next_id += 1
    return NonMarkovCell(pop._next_id, t, j, f)
end

"""
    celldeath!(population, node)

Remove the dead cell `node` from the population and prune it from the lineage tree.
"""
function celldeath!(
    population::Population,
    node::BinaryNode{NonMarkovCell},
)
    prune_tree!(node)
    delete!(population.cells, node.data.id)
end
