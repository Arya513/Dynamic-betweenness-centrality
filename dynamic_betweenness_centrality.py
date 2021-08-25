import networkx as nx


# TODO: check restrictions s != v != t when calculating A, B, C, A', B' and C'

# version of the algorithm where only insertion of edges is allowed
def icentral_incremental(G: nx.Graph, from_: int, to_: int, bc: dict) -> dict:
    # add egde -> make new graph
    G_with_added_edge = G.copy()
    G_with_added_edge.add_edge(from_, to_)

    # get list of biconnected components in a list
    # this has to be done after edge insertion (in case two biconnected components merge into one)
    biconnected_components = list(nx.biconnected_components(G_with_added_edge))

    # these edges are needed to determine which biconnected component is affected
    # in other words - to find to which biconnected component the inserted edge belongs to
    biconnected_components_edges = list(nx.biconnected_component_edges(G_with_added_edge))

    # find index of affected component
    index_of_affected_biconnected_component: int
    for i, edges_of_component in enumerate(biconnected_components_edges):
        if (from_, to_) in edges_of_component or (to_, from_) in edges_of_component:
            index_of_affected_biconnected_component = i

    # get subgraph of affected biconnected component without added edge
    Be_subgraph_without_edge = G.subgraph(biconnected_components[index_of_affected_biconnected_component])

    # get subgraph of affected biconnected component with added edge
    Be_subgraph_with_edge = G_with_added_edge.subgraph(biconnected_components[index_of_affected_biconnected_component])

    # TODO: not sure how this part would work in directed graphs
    # we make two BFS' within affected component - one starting from from_ node, the other starting from to_
    # we keep the distances from the source node in the two dictionaries below
    bfs1 = bfs_counting_hops(Be_subgraph_without_edge, from_)
    bfs2 = bfs_counting_hops(Be_subgraph_without_edge, to_)

    # we create the set q
    q = set()
    for node in bfs1.keys():
        if bfs1[node] != bfs2[node]:
            q.add(node)

    # get articulation points
    articulation_points = list(nx.algorithms.components.articulation_points(G))

    # PART 1 (lines 11 - 25 of pseudocode)
    for s in q:
        # create sigma_s and predecessors_s
        sigma_s = dict()
        predecessors_s = dict()

        # sigma_s - the starting node is s and the target nodes (t) are all nodes IN Be
        # predecessors_s - predecessors of a node on the shortest path between s (source) and current node (t)
        # we want to count how many shortest paths are between s and t
        for t in Be_subgraph_without_edge.nodes:
            sigma_s[t] = len(list(nx.algorithms.shortest_paths.generic.all_shortest_paths(G, source=s, target=t)))

            # TODO: extend for directed graphs
            # TODO: check if s is included as predecessor
            predecessors_s[t] = nx.algorithms.shortest_paths.unweighted.predecessor(G, source=s, target=t)

        # initialize delta_s(v) and delta_Gs(v)
        delta_s = dict()
        delta_Gs = dict()

        for v in Be_subgraph_without_edge.nodes:
            delta_Gs[v] = 0
            delta_s[v] = 0

        # get nodes in reverse BFS order from source node s
        nodes_in_reverse_bfs_order = get_reverse_BFS_order(Be_subgraph_without_edge, s)

        # TODO: check if reverse BFS order excludes s
        # nodes_in_reverse_bfs_order.remove(s)

        for w in nodes_in_reverse_bfs_order:
            if s in articulation_points and w in articulation_points:
                delta_Gs[w] = get_cardinality_of_Gi(G, s,
                                                    biconnected_components[index_of_affected_biconnected_component]) * \
                              get_cardinality_of_Gi(G, w,
                                                    biconnected_components[index_of_affected_biconnected_component])

            for p in predecessors_s[w]:
                delta_s[p] = delta_s[p] + float(sigma_s[p]) / float(sigma_s[w]) * (1 + delta_s[w])

                if s in articulation_points:
                    delta_Gs[p] = delta_Gs[p] + delta_Gs[w] * float(sigma_s[p]) / float(sigma_s[w])

            if w != s:
                bc[w] = bc[w] - float(delta_s[w]) / 2.0

            if s in articulation_points:
                bc[w] = bc[w] - delta_s[w] * get_cardinality_of_Gi(G, s,
                                                                   biconnected_components[
                                                                       index_of_affected_biconnected_component])
                bc[w] = bc[w] - float(delta_Gs[w]) / 2.0

        # PART 2 (lines 26 - 40 of pseudocode)

        # create sigma_s2 and predecessors_s2 (exactly the same as in the first part of the algorithm - the only
        # difference is the added edge)
        sigma_s2 = dict()
        predecessors_s2 = dict()

        # sigma_s2 - the starting node is s and the target nodes (t) are all nodes IN Be' (Be' == Be for vertices)
        # predecessors_s2 - predecessors of a node on shortest path between s (source) and current node (t)
        # we want to count how many shortest paths are between s and t
        for t in Be_subgraph_with_edge.nodes:
            sigma_s2[t] = len(list(nx.algorithms.shortest_paths.generic.all_shortest_paths(G_with_added_edge,
                                                                                           source=s,
                                                                                           target=t)))

            # TODO: extend for directed graphs
            # TODO: check if s is included as predecessors
            predecessors_s2[t] = nx.algorithms.shortest_paths.unweighted.predecessor(G_with_added_edge,
                                                                                     source=s,
                                                                                     target=t)

        # initialize delta_s2(v) and delta_Gs2(v)
        delta_s2 = dict()
        delta_Gs2 = dict()

        for v in Be_subgraph_with_edge.nodes:
            delta_Gs2[v] = 0
            delta_s2[v] = 0

        # get nodes in reverse BFS order from source node s
        nodes_in_reverse_bfs_order = get_reverse_BFS_order(Be_subgraph_with_edge, s)

        # TODO: check if reverse BFS order excludes s
        # nodes_in_reverse_bfs_order.remove(s)

        for w in nodes_in_reverse_bfs_order:
            if s in articulation_points and w in articulation_points:
                delta_Gs2[w] = get_cardinality_of_Gi(G_with_added_edge, s,
                                                     biconnected_components[
                                                         index_of_affected_biconnected_component]) * \
                               get_cardinality_of_Gi(G_with_added_edge, w,
                                                     biconnected_components[
                                                         index_of_affected_biconnected_component])

            for p in predecessors_s2[w]:
                delta_s2[p] = delta_s2[p] + float(sigma_s2[p]) / float(sigma_s2[w]) * (1 + delta_s2[w])

                if s in articulation_points:
                    delta_Gs2[p] = delta_Gs2[p] + delta_Gs2[w] * float(sigma_s2[p]) / float(sigma_s2[w])

            if w != s:
                bc[w] = bc[w] + float(delta_s2[w]) / 2.0

            if s in articulation_points:
                bc[w] = bc[w] + delta_s2[w] * get_cardinality_of_Gi(G_with_added_edge, s,
                                                                    biconnected_components[
                                                                        index_of_affected_biconnected_component])
                bc[w] = bc[w] + float(delta_Gs2[w]) / 2.0

    return bc


# get number of nodes in subgraph connected to Be via articulation point a
def get_cardinality_of_Gi(G: nx.Graph, a: int, V_Be: list) -> int:
    visited = list()
    q = list()

    count = 0

    visited.append(a)
    q.append(a)

    while q:
        a = q.pop(0)
        for n in G.neighbors(a):
            if n not in V_Be:
                if n not in visited:
                    visited.append(n)
                    q.append(n)
                    count += 1

    return count


# function for getting reverse BFS order
def get_reverse_BFS_order(subgraph: nx.Graph, s: int) -> list:
    visited = list()
    q = list()

    visited.append(s)
    q.append(s)

    while q:
        s = q.pop(0)
        for n in subgraph.neighbors(s):
            if n not in visited:
                visited.append(n)
                q.append(n)

    visited.reverse()

    return visited


# version of BFS counting and returning number of hops from node to source node
def bfs_counting_hops(G: nx.Graph, s: int) -> dict:
    distances = dict()

    distances[s] = 0

    visited = list()
    q = list()

    visited.append(s)
    q.append((s, 0))

    while q:
        s, distance = q.pop(0)

        for n in G.neighbors(s):

            if n not in visited:
                visited.append(n)
                q.append((n, distance + 1))
                distances[n] = distance + 1

    return distances


# example from paper
def make_graph() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(1, 21))
    G.add_edge(9, 10)
    G.add_edge(9, 11)
    G.add_edge(10, 1)
    G.add_edge(11, 1)
    G.add_edge(1, 2)
    G.add_edge(1, 4)
    G.add_edge(2, 3)
    G.add_edge(4, 3)
    G.add_edge(4, 6)
    G.add_edge(6, 7)
    G.add_edge(7, 8)
    G.add_edge(6, 12)
    G.add_edge(6, 14)
    G.add_edge(12, 13)
    G.add_edge(14, 13)
    G.add_edge(3, 5)
    G.add_edge(8, 5)
    G.add_edge(5, 15)
    G.add_edge(5, 19)
    G.add_edge(19, 17)
    G.add_edge(15, 17)
    G.add_edge(17, 16)
    G.add_edge(17, 20)
    G.add_edge(16, 18)
    G.add_edge(20, 18)

    return G


def main():
    G = make_graph()

    G_added_edge = make_graph()
    G_added_edge.add_edge(4, 8)

    bc_brandes = nx.betweenness_centrality(G_added_edge)
    bc_incentral = icentral_incremental(G, 4, 8, nx.betweenness_centrality(G))

    # for node in bc_incentral.keys():
    #    if bc_brandes[node] != bc_incentral[node]:
    #        print("Failed")
    #        break

    for k, v in bc_incentral.items():
        print(f"{k}: {v}")


if __name__ == '__main__':
    main()
