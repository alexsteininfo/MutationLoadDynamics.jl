"""
    simulate!(population, block, rng) -> Population

Run a non-Markovian birth-death simulation on `population` using `block` parameters.
Each living cell pre-schedules its next event (division or death) at birth by drawing
from competing gamma (or arbitrary) distributions — the earlier of the two fires.
A global min-heap orders events so the next event is always processed first.

Returns `population` (mutated in-place) so blocks can be chained.
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

    while true
        heap = BinaryMinHeap{CellEvent}()
        for node in values(population.cells)
            schedule_cell!(heap, node, block, rng)
        end

        isnothing(accumulator) || record_trajectory_if_due!(accumulator, population)

        while !block.stopfunction(population) && !isempty(heap)
            event = pop!(heap)
            population.t = event.time

            if event.event_type == :birth
                d1, d2 = celldivision!(population, event.node, event.time, block, rng)
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
            isnothing(accumulator) || _reset_accumulator!(accumulator)
        else
            break
        end
    end

    isnothing(accumulator) || _fire_end_triggers!(accumulator, population)
    return population
end
