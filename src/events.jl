# `CellEvent` and its `Base.isless` method are defined in `types.jl`, because
# `Population` needs the type in order to hold a queue of pending events.

# Rejection draws allowed before `schedule_cell!` gives up on conditioning an old cell.
const _MAX_CONDITIONING_DRAWS = 100_000

"""
    schedule_cell!(heap, node, block, rng, tmin = node.data.birthtime)

Pre-schedule the next event for `node` using competing risks: draw a division time from
`block.birth_dist` and a death time from `block.death_dist`, both measured from the
cell's birth, and push the earlier one onto `heap`.

`tmin` is the earliest time the event may fire. For a newborn cell it is the birthtime
and the first draw is always accepted. When an *existing* cell is rescheduled (a fresh
population with old birthtimes, [`reset_schedule!`](@ref), or an extinction restart),
`tmin` is the current time and both draws are repeated until the earlier one falls at or
after `tmin`. That is exact rejection sampling from the law of the cell's next event
given that nothing happened to it before `tmin`.
"""
function schedule_cell!(
    heap::BinaryMinHeap{CellEvent},
    node::BinaryNode{NonMarkovCell},
    block::NonMarkovBlock,
    rng::AbstractRNG,
    tmin::Float64 = node.data.birthtime,
)
    f  = node.data.fitness
    t0 = node.data.birthtime
    for _ in 1:_MAX_CONDITIONING_DRAWS
        t_div = t0 + rand(rng, block.birth_dist(f))
        t_die = t0 + rand(rng, block.death_dist(f))
        if min(t_div, t_die) >= tmin
            if t_div <= t_die
                push!(heap, CellEvent(t_div, node, :birth))
            else
                push!(heap, CellEvent(t_die, node, :death))
            end
            return nothing
        end
    end
    error("cell $(node.data.id), born at t = $t0, could not be given an event at or " *
          "after t = $tmin in $_MAX_CONDITIONING_DRAWS draws: under this block it would " *
          "almost surely have divided or died already. Its age is incompatible with the " *
          "waiting-time distributions.")
end
