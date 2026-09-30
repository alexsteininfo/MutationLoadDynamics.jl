"""
    prune_tree!(cellnode)

Remove `cellnode` from the tree and recursively remove any ancestor that becomes
childless. Called when a leaf cell dies without dividing.
"""
function prune_tree!(cellnode)
    while true
        parent = cellnode.parent
        if isnothing(parent)
            return
        else
            cellnode.parent = nothing
            if parent.left === cellnode
                parent.left = nothing
            elseif parent.right === cellnode
                parent.right = nothing
            else
                error("dead cell is neither left nor right child of parent")
            end
            if isnothing(parent.left) && isnothing(parent.right)
                cellnode = parent
            else
                return
            end
        end
    end
end

"""
    endtime(cellnode::BinaryNode)

Return the time at which the cell divided (birthtime of its left child), or `nothing`
if the cell is still alive (a leaf).
"""
function endtime(cellnode::BinaryNode)
    if !isnothing(cellnode.left)
        return cellnode.left.data.birthtime
    elseif !isnothing(cellnode.right)
        return cellnode.right.data.birthtime
    else
        return nothing
    end
end

"""
    celllifetime(cellnode::BinaryNode, [tmax])

Compute the lifetime of a cell. If it has not yet divided, use `tmax` as the end time
(defaults to the age of the tree).
"""
function celllifetime(cellnode::BinaryNode, tmax=nothing)
    et = endtime(cellnode)
    if !isnothing(et)
        return et - cellnode.data.birthtime
    else
        tmax = isnothing(tmax) ? age(getroot(cellnode)) : tmax
        return tmax - cellnode.data.birthtime
    end
end

"""
    celllifetimes(root; excludeliving=true)

Compute the lifetime of every cell in the phylogeny rooted at `root`, in pre-order.
By default, currently alive (leaf) cells are excluded.
"""
function celllifetimes(root::BinaryNode; excludeliving::Bool = true)
    nodes, _ = _preorder(root)
    if excludeliving
        return Float64[endtime(n) - n.data.birthtime for n in nodes if haschildren(n)]
    else
        popage = age(root)
        return Float64[celllifetime(n, popage) for n in nodes]
    end
end

"""
    age(root::BinaryNode)

Return the birthtime of the most recently born leaf cell, i.e. the time of the last
division in the tree. This is not `pop.t` if the run ended on a death.
"""
age(root::BinaryNode) = maximum(leaf.data.birthtime for leaf in _leaves(root))

"""
    age(population::Population)

Return `population.t`.
"""
age(population::Population) = population.t

"""
    getalivecells(root::BinaryNode) -> Vector

Return all alive leaf cells descending from `root`, in `Leaves(root)` order. Dead cells
are pruned from the tree, so every leaf is alive.
"""
getalivecells(root::BinaryNode) = _leaves(root)
