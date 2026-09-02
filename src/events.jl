# `CellEvent` and its `Base.isless` method are defined in `types.jl`, because
# `Population` needs the type in order to hold a queue of pending events.

"""
    schedule_cell!(heap, node, block, rng)

Pre-schedule the next event for `node` using the competing-risks approach:
draw `t_div` from `block.birth_dist` and `t_die` from `block.death_dist`,
then push the earlier event onto `heap`.
"""
function schedule_cell!(
    heap::BinaryMinHeap{CellEvent},
    node::BinaryNode{NonMarkovCell},
    block::NonMarkovBlock,
    rng::AbstractRNG,
)
    f  = node.data.fitness
    t0 = node.data.birthtime
    t_div = t0 + rand(rng, block.birth_dist(f))
    t_die = t0 + rand(rng, block.death_dist(f))
    if t_div <= t_die
        push!(heap, CellEvent(t_div, node, :birth))
    else
        push!(heap, CellEvent(t_die, node, :death))
    end
end
