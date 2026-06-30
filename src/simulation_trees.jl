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
    if haschildren(cellnode)
        return cellnode.left.data.birthtime
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
    if haschildren(cellnode)
        return cellnode.left.data.birthtime - cellnode.data.birthtime
    else
        tmax = isnothing(tmax) ? age(getroot(cellnode)) : tmax
        return tmax - cellnode.data.birthtime
    end
end

"""
    celllifetimes(root; excludeliving=true)

Compute the lifetime of every cell in the phylogeny rooted at `root`.
By default, currently alive (leaf) cells are excluded.
"""
function celllifetimes(root; excludeliving=true)
    lifetimes = Float64[]
    if excludeliving
        for cellnode in PreOrderDFS(root)
            if haschildren(cellnode)
                push!(lifetimes, cellnode.left.data.birthtime - cellnode.data.birthtime)
            end
        end
    else
        popage = age(root)
        for cellnode in PreOrderDFS(root)
            push!(lifetimes, celllifetime(cellnode, popage))
        end
    end
    return lifetimes
end

"""
    age(root::BinaryNode)

Return the birthtime of the most recently born leaf cell — the current simulation age.
"""
function age(root::BinaryNode)
    t = 0.0
    for cellnode in Leaves(root)
        if cellnode.data.birthtime > t
            t = cellnode.data.birthtime
        end
    end
    return t
end

"""
    age(population::Population)

Return `population.t`.
"""
age(population::Population) = population.t

isalive(cellnode::BinaryNode{T}) where T = isalive(cellnode.data)
isalive(cell::NonMarkovCell) = true
isalive(::Nothing) = false

"""
    getalivecells(root::BinaryNode) -> Vector

Return all alive leaf cells descending from `root`.
"""
getalivecells(root::BinaryNode) =
    [cellnode for cellnode in Leaves(root) if isalive(cellnode.data)]

"""
    asroot!(node)

Temporarily detach `node` from its parent by setting `parent = nothing`.
Returns `(node, original_parent)`.
"""
function asroot!(node)
    parent = node.parent
    node.parent = nothing
    return node, parent
end
