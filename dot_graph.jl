module DotGraph
using GraphViz

struct Node
    id :: Int
    label :: String
end

struct Edge
    from :: Int
    to :: Int
    label :: String
end

mutable struct Graph
    directed :: Bool
    fresh_node_idx :: Int
    nodes :: Vector{Node}
    edges :: Vector{Edge}

    Graph(directed :: Bool) = new(directed, 1, [], [])
end

function add_node!(graph :: Graph, label :: String) :: Int 
    id = graph.fresh_node_idx
    graph.fresh_node_idx += 1
    push!(graph.nodes, Node(id, label))
    return id
end

function add_edge!(graph :: Graph, from :: Int, to :: Int, label :: String = "")
    push!(graph.edges, Edge(from, to, label))
    return
end

# Build the DOT string corresponding to the graph.
function build(graph :: Graph) :: String
    buf = IOBuffer()
    if graph.directed
        write(buf, "digraph {\n")
    else
        write(buf, "graph {\n")
    end

    # Add nodes.
    for node in graph.nodes
        write(buf, "    $(node.id) [ label=\"$(node.label)\"]\n")
    end

    # Add edges.
    for edge in graph.edges
        arrow = graph.directed ? "->" : "--"
        if edge.label == ""
            write(buf, "    $(edge.from) $(arrow) $(edge.to)\n")
        else 
            write(buf, "    $(edge.from) $(arrow) $(edge.to) [ label=\"$(edge.label)\" ]\n")
        end
    end

    write(buf, "}\n")
    return String(take!(buf))
end

function render(graph :: Graph)
    # For some reason I have to rebuild the buffer here,
    # it crashes if I use the buffer of [build] directly.
    return GraphViz.load(IOBuffer(build(graph)))
end

end
