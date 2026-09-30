"""
    simulate!(population, block, rng; accumulator = nothing) -> Population

Run a non-Markovian birth-death simulation on `population` using `block` parameters.
Each living cell pre-schedules its next event (division or death) at birth by drawing
from competing waiting-time distributions — the earlier of the two fires. A global
min-heap orders events so the globally earliest one is always processed next.

Returns `population` (mutated in place) so blocks can be chained.

# Chaining

Calling `simulate!` again on the same `population` continues the same lineage tree.
The queue of already-drawn, not-yet-fired events is carried on `population._pending`
and reused, so a chained run is identical draw for draw to the equivalent uninterrupted
one. A carried event was drawn under the *previous* block's `birth_dist`/`death_dist`;
call [`reset_schedule!`](@ref) first to redraw every cell under the new block instead.

# Fresh schedules

Whenever there is no carried queue — a new population, after `reset_schedule!`, or
after an extinction restart — every alive cell is scheduled in id order, conditioned
on its current age: its next event is drawn from the law of `(T_div, T_die)` given that
neither happened before `population.t`
(see [`schedule_cell!`](@ref MutationLoadDynamics.schedule_cell!)). This is exact, so
restarts and redrawn schedules are also correct on chained, aged populations.

The rng defaults to `Random.default_rng()`; pass your own for reproducible runs.
"""
function simulate!(
    population::Population,
    block::NonMarkovBlock,
    rng::AbstractRNG = Random.default_rng();
    accumulator::Union{MeasurementAccumulator, Nothing} = nothing,
)
    if block.restart_on_extinction
        # A deepcopy reaches the whole tree through `parent` links, so this is cheap only
        # while the population is small — typically on a first block.
        initial_cells  = deepcopy(population.cells)
        initial_t      = population.t
        initial_nextid = population._next_id
    end
    drivers = Poisson(block.ν)

    if !isnothing(accumulator)
        _start_call!(accumulator, population)
        checkpoint = _checkpoint(accumulator)
    end

    while true
        heap = population._pending
        if isnothing(heap)
            heap = BinaryMinHeap{CellEvent}()
            for id in sort!(collect(keys(population.cells)))
                schedule_cell!(heap, population.cells[id], block, rng, population.t)
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

        isnothing(accumulator) || _check_popsize_triggers!(accumulator, population)

        while !block.stopfunction(population) && !isempty(heap)
            event = pop!(heap)
            event.time >= population.t || error(
                "event at t = $(event.time) precedes the current time t = " *
                "$(population.t); the event queue is corrupt")
            # Record grid points and AtTime snapshots due before this event, while the
            # state is still exactly the state that held at those times.
            isnothing(accumulator) || _record_until!(accumulator, population, event.time)
            population.t = event.time

            if event.event_type == :birth
                d1, d2 = celldivision!(population, event.node, event.time, block,
                                       drivers, rng)
                # Must run before scheduling: schedule_cell! reads the fitness at push
                # time, so a boosted daughter's own first division uses the new fitness.
                isnothing(block.on_division) ||
                    block.on_division(population, event.node, d1, d2)
                schedule_cell!(heap, d1, block, rng)
                schedule_cell!(heap, d2, block, rng)
            else
                celldeath!(population, event.node)
            end

            isnothing(accumulator) || _check_popsize_triggers!(accumulator, population)
        end

        if popsize(population) == 0 && block.restart_on_extinction
            population.cells    = deepcopy(initial_cells)
            population.t        = initial_t
            population._next_id = initial_nextid
            # deepcopy produces fresh BinaryNodes, so any carried events point at
            # orphaned nodes — drop the queue and redraw for the restored cells.
            population._pending = nothing
            isnothing(accumulator) || _rollback!(accumulator, checkpoint)
            isnothing(block.on_restart) || block.on_restart(population)
        else
            break
        end
    end

    isnothing(accumulator) || _finish_call!(accumulator, population)
    return population
end

"""
    reset_schedule!(population) -> Population

Discard the queue of pending, already-drawn events carried on `population`, so that the
next `simulate!` call redraws every alive cell's next event from scratch.

This is the escape hatch from the event-carrying behaviour described under [`simulate!`](@ref).
It is also how to recover after adding or removing cells in `population.cells` by hand,
which otherwise leaves the queue inconsistent with the population and raises an error.

The redraw conditions each cell on the age it has already reached (see
[`simulate!`](@ref)), so it is exact; the only change is that every cell's next event
now follows the distributions of the block passed to the next `simulate!` call.
"""
function reset_schedule!(population::Population)
    population._pending = nothing
    return population
end
