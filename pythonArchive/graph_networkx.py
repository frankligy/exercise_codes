import networkx as nx

G = nx.Graph()  # networkx.classes.graph.Graph
G = nx.DiGraph()  # networkx.classes.digraph.DiGraph
G = nx.MultiGraph()  # networkx.classes.multigraph.MultiGraph
G = nx.MultiDiGraph()  # networkx.classes.multidigraph.MultiDiGraph

# Graph
G = nx.Graph()
G.add_nodes_from([1,2,3,4,5],color='blue')
G.add_edges_from([(1,2,{'weight':2}),(4,5,{'weight':5})])

G.nodes(data=True)
G.nodes.data('color')
G.nodes[1]
G.edges(data=True)
G.edges.data('weight')
G.edges[1,2]
G.adj

#



