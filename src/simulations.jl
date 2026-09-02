"""
    simulate!(population, block, rng) -> Population

Run a non-Markovian birth-death simulation on `population` using `block` parameters.
Each living cell pre-schedules its next event (division or death) at birth by drawing
from competing gamma (or arbitrary) distributions — the earlier of the two fires.
A global min-heap orders events so the next event is always processed first.

Returns `population` (mutated in-place) so blocks can be chained.

# Chaining

Calling `simulate!` again on the same `population` continues the same lineage tree.
The queue of already-drawn, not-yet-fired events is carried on `population._pending`
and **reused**, rather than redrawn.

This matters for any non-exponential waiting-time distribution. Redrawing would anchor
each cell's next event at its `birthtime`, ignoring the age it has already accumulated
without an event; the correct law is the conditional `T | T > age`, and drawing
unconditionally both loses that conditioning and can place events before the current
`population.t`. Reusing the pending events sidesteps the problem exactly.

The consequence to be aware of: a carried event was drawn under the *previous* block's
`birth_dist`/`death_dist`. The cell has already committed to that division or death time,
so a second block that changes those distributions only affects cells scheduled after the
boundary. Call [`reset_schedule!`](@ref) first if you want every cell redrawn under the
new law instead — but note that this reintroduces the age-conditioning error above, so it
is only sound when all cells were just born.

Chaining is not supported together with `restart_on_extinction`: the restored cells carry
birthtimes from before the restart point, so their rescheduling has the same defect.
"""
function simulate!(
    population::Population,
    block::NonMarkovBlock,
    rng::AbstractRNG = Random.GLOBAL_RNG;
    accumulator::Union{MeasurementAccumulator, Nothing} = nothing,
)
    if block.restart_on_extinction
        initial_cells  = deepcopy(population.cells)
        initial_t      = population.t
        initial_nextid = population._next_id
    end

    # Start trajectory recording at the current time. Without this a fresh accumulator
    # handed to a chained call would back-fill points from t = 0 up to population.t.
    if !isnothing(accumulator)
        accumulator.next_trajectory_t =
            max(accumulator.next_trajectory_t, population.t)
    end

    while true
        heap = population._pending
        if isnothing(heap)
            heap = BinaryMinHeap{CellEvent}()
            for node in values(population.cells)
                schedule_cell!(heap, node, block, rng)
            end
        elseif length(heap) != popsize(population)
            throw(ArgumentError(
                "population's carried event queue holds $(length(heap)) events but the " *
                "population has $(popsize(population)) alive cells — cells were added or " *
                "removed outside simulate!. Call reset_schedule!(population) to discard " *
                "the queue and redraw every cell's next event from scratch."))
        end
        # BinaryMinHeap mutates in place, so this reference stays current as events fire.
        population._pending = heap

        isnothing(accumulator) || record_trajectory_if_due!(accumulator, population)

        while !block.stopfunction(population) && !isempty(heap)
            event = pop!(heap)
            population.t = event.time

            if event.event_type == :birth
                d1, d2 = celldivision!(population, event.node, event.time, block, rng)
                # Must run before scheduling: schedule_cell! reads node.data.fitness at
                # push time, so a hook placed after it would leave a boosted daughter's
                # own first division drawn at the pre-boost fitness.
                isnothing(block.on_division) ||
                    block.on_division(population, event.node, d1, d2)
                schedule_cell!(heap, d1, block, rng)
                schedule_cell!(heap, d2, block, rng)
            else
                celldeath!(population, event.node)
            end

            if !isnothing(accumulator)
                record_trajectory_if_due!(accumulator, population)
                check_timed_triggers!(accumulator, population)
            end
        end

        N = popsize(population)
        if N == 0 && block.restart_on_extinction
            population.cells    = deepcopy(initial_cells)
            population.t        = initial_t
            population._next_id = initial_nextid
            # deepcopy produces fresh BinaryNodes, so any carried events point at
            # orphaned nodes — drop the queue and redraw for the restored cells.
            population._pending = nothing
            isnothing(accumulator) || _reset_accumulator!(accumulator, initial_t)
            isnothing(block.on_restart) || block.on_restart(population)
        else
            break
        end
    end

    isnothing(accumulator) || _fire_end_triggers!(accumulator, population)
    return population
end

"""
    reset_schedule!(population) -> Population

Discard the queue of pending, already-drawn events carried on `population`, so that the
next `simulate!` call redraws every alive cell's next event from scratch.

This is the escape hatch from the event-carrying behaviour described under [`simulate!`](@ref).
It is also how to recover after adding or removing cells in `population.cells` by hand,
which otherwise leaves the queue inconsistent with the population and raises an error.

Be aware that redrawing anchors each cell at its `birthtime` and ignores the age it has
already accumulated, so it is only statistically sound when every alive cell has just been
born (as in a freshly initialised population).
"""
function reset_schedule!(population::Population)
    population._pending = nothing
    return population
end
