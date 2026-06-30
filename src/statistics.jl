"""
    mutations_per_cell(population::Population) -> Vector{Int64}

Total accumulated driver mutations per alive cell, summing each cell's lineage back
to the root.
"""
function mutations_per_cell(population::Population)
    mutspercell = Int64[]
    for cellnode in allcells(population)
        mutations = cellnode.data.mutations
        node = cellnode
        while !isnothing(node.parent)
            node = node.parent
            mutations += node.data.mutations
        end
        push!(mutspercell, mutations)
    end
    return mutspercell
end

function mutations_per_cell(root::BinaryNode{T}; includeclonal=false) where T <: AbstractTreeCell
    mutspercell = Int64[]
    for cellnode in Leaves(root)
        if isalive(cellnode.data)
            mutations = cellnode.data.mutations
            node = cellnode
            while true
                if !AbstractTrees.isroot(node) && (node != root || includeclonal)
                    node = node.parent
                    mutations += node.data.mutations
                else
                    break
                end
            end
            push!(mutspercell, mutations)
        end
    end
    return mutspercell
end

"""
    average_mutations(population::Population) -> Float64

Mean accumulated driver mutations across all alive cells.
"""
average_mutations(population::Population) = mean(mutations_per_cell(population))

"""
    clonal_mutations(population::Population) -> Int64

Number of driver mutations shared by every alive cell (mutations at the MRCA).
"""
function clonal_mutations(population::Population)
    MRCA = findMRCA(population)
    isnothing(MRCA) && return 0
    return all_cell_mutations(MRCA)
end

"""
    all_cell_mutations(node::BinaryNode) -> Int64

Total mutations at `node` including all inherited ones (sum to root).
"""
function all_cell_mutations(node::BinaryNode)
    mutations = node.data.mutations
    while !isnothing(node.parent)
        node = node.parent
        mutations += node.data.mutations
    end
    return mutations
end

"""
    fitness_per_cell(population::Population) -> Vector{Float64}

Fitness of each currently alive cell (stored directly in `NonMarkovCell.fitness`).
"""
fitness_per_cell(population::Population) =
    [node.data.fitness for node in values(population.cells)]

"""
    fitness_distribution(population::Population) -> Vector{Float64}

Alias for `fitness_per_cell`; returns fitness values of all living cells.
"""
fitness_distribution(population::Population) = fitness_per_cell(population)

"""
    mean_k(population::Population) -> Float64

Mean total driver-mutation count per alive cell.
"""
mean_k(population::Population) = mean(Float64.(mutations_per_cell(population)))

"""
    var_k(population::Population) -> Float64

Variance in total driver-mutation count across alive cells.
"""
var_k(population::Population) = var(Float64.(mutations_per_cell(population)))

"""
    pairwise_differences(population::Population[, idx]) -> Dict{Int64,Int64}

Histogram of pairwise mutational distances between alive cells.
"""
function pairwise_differences(population::Population, idx=nothing)
    cells = allcells(population)
    isnothing(idx) || (cells = cells[idx])
    n = length(cells)
    pfd_vec = Int64[]
    for i in 1:n, j in i+1:n
        push!(pfd_vec, pairwisedistance(cells[i], cells[j]))
    end
    return countmap(pfd_vec)
end

"""
    pairwisedistance(node1, node2) -> Int64

Number of driver mutations that differ between two cells (walk to their MRCA).
"""
function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode)
    cellnode1 == cellnode2 && return 0
    distance = 0
    while true
        if cellnode1.data.id > cellnode2.data.id
            cellnode1, cellnode2 = cellnode2, cellnode1
        end
        if (cellnode1.parent == cellnode2.parent) ||
           (isnothing(cellnode1.parent) && isnothing(cellnode2.parent))
            return distance + cellnode1.data.mutations + cellnode2.data.mutations
        elseif cellnode1 == cellnode2.parent
            return distance + cellnode2.data.mutations
        else
            distance += cellnode2.data.mutations
            cellnode2 = cellnode2.parent
        end
    end
end

"""
    pairwisedistances(population::Population[, idx]) -> Vector{Int64}

All pairwise distances between alive cells as a flat vector.
"""
function pairwisedistances(population::Population, idx=nothing)
    cells = allcells(population)
    isnothing(idx) || (cells = cells[idx])
    n = length(cells)
    pfd_vec = Int64[]
    for i in 1:n, j in i+1:n
        push!(pfd_vec, pairwisedistance(cells[i], cells[j]))
    end
    return pfd_vec
end

function _time_to_MRCA(node1, node2, t)
    if node1.data.birthtime > node2.data.birthtime
        node1, node2 = node2, node1
    end
    isdefined(node1, :parent) || return t - endtime(node1)
    isdefined(node2, :parent) || return t - endtime(node2)
    if node1.parent == node2.parent
        return t - node1.data.birthtime
    else
        return _time_to_MRCA(node1, node2.parent, t)
    end
end

"""
    coalescence_times(root[, idx]; t) -> Vector{Float64}
    coalescence_times(population[, idx]; t) -> Vector{Float64}

Time to MRCA (coalescence time) for every pair of alive cells.
"""
function coalescence_times(root::BinaryNode, idx=nothing; t=nothing)
    t = isnothing(t) ? age(root) : t
    coaltimes = Float64[]
    alivecells = getalivecells(root)
    isnothing(idx) || (alivecells = alivecells[idx])
    while length(alivecells) > 1
        c1 = popfirst!(alivecells)
        for c2 in alivecells
            push!(coaltimes, _time_to_MRCA(c1, c2, t))
        end
    end
    return coaltimes
end

function coalescence_times(population::Population, idx=nothing; t=nothing)
    t = isnothing(t) ? age(population) : t
    coaltimes = Float64[]
    alivecells = allcells(population)
    isnothing(idx) || (alivecells = alivecells[idx])
    while length(alivecells) > 1
        c1 = popfirst!(alivecells)
        for c2 in alivecells
            push!(coaltimes, _time_to_MRCA(c1, c2, t))
        end
    end
    return coaltimes
end

"""
    sitefrequencyspectrum(population::Population) -> Vector{Int64}

Driver-mutation site-frequency spectrum: `sfs[k]` = number of driver mutation events
present in exactly `k` living cells.
"""
function sitefrequencyspectrum(population::Population)
    sfs  = zeros(Int64, popsize(population))
    root = getsingleroot(allcells(population))
    if isnothing(root)
        # multiple independent roots (forest) — traverse each root separately
        roots = unique(n -> objectid(_treeroot(n)), values(population.cells))
        for r in roots
            _sfs_fill!(r, sfs)
        end
    else
        _sfs_fill!(root, sfs)
    end
    return sfs
end

function _treeroot(node::BinaryNode)
    while !isnothing(node.parent)
        node = node.parent
    end
    return node
end

function _sfs_fill!(node::BinaryNode, sfs::Vector{Int64})
    if isnothing(node.left) && isnothing(node.right)
        sfs[1] += node.data.mutations
        return 1
    end
    count = 0
    isnothing(node.left)  || (count += _sfs_fill!(node.left,  sfs))
    isnothing(node.right) || (count += _sfs_fill!(node.right, sfs))
    count > 0 && (sfs[count] += node.data.mutations)
    return count
end
