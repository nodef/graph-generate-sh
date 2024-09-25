import numpy as np
import scipy.sparse as sp
import os
import matplotlib.pyplot as plt
from scipy.sparse.linalg import svds
import sys
from collections import deque

def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

def read_graph(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    n, m = map(int, lines[0].split())
    
    adj_matrix = sp.lil_matrix((n, n), dtype=np.float64)
    
    for line in lines[1:]:
        node1, node2, wt = map(int, line.split())
        adj_matrix[node1-1, node2-1] = wt+1  
    
    return adj_matrix, n


def find_bridges(adj_matrix,component, n):
    """Finds bridges in the graph."""
    bridges = []
    visited = [False] * n
    discovery = [-1] * n
    low = [-1] * n
    parent = [-1] * n
    time = [0]  # This needs to be a list so it's mutable across recursion
    
    def bridge_util(u):
        visited[u] = True
        discovery[u] = low[u] = time[0]
        time[0] += 1
        
        for v in adj_matrix.rows[u]:
            if not visited[v]:
                parent[v] = u
                bridge_util(v)
                
                low[u] = min(low[u], low[v])
                
                # If the lowest vertex reachable from v is below u in DFS tree,
                # then u-v is a bridge
                if low[v] > discovery[u]:
                    bridges.append((u, v))
            elif v != parent[u]:
                low[u] = min(low[u], discovery[v])
    
    for i in component:
        if not visited[i]:
            bridge_util(i)
    
    return bridges

def find_sccs(adj_matrix, n):
    """Find all SCCs using Kosaraju's algorithm."""
    visited = [False] * n
    stack = []
    
    def fill_order(v):
        visited[v] = True
        for u in adj_matrix.rows[v]:
            if not visited[u]:
                fill_order(u)
        stack.append(v)
    
    def dfs(v, transposed_adj, visited, component):
        visited[v] = True
        component.append(v)
        for u in transposed_adj.rows[v]:
            if not visited[u]:
                dfs(u, transposed_adj, visited, component)
    
    for i in range(n):
        if not visited[i]:
            fill_order(i)
    
    transposed_adj = adj_matrix.transpose()
    
    visited = [False] * n
    sccs = []
    while stack:
        v = stack.pop()
        if not visited[v]:
            component = []
            dfs(v, transposed_adj, visited, component)
            sccs.append(component)
    
    return sccs


def bfs(capacity, source, sink, parent):
    visited = [False] * len(capacity)
    queue = deque([source])
    visited[source] = True

    while queue:
        u = queue.popleft()
        for v in range(len(capacity)):
            if not visited[v] and capacity[u][v] > 0:
                queue.append(v)
                visited[v] = True
                parent[v] = u
                if v == sink:
                    return True
    return False

def edmonds_karp(capacity, source, sink):
    """Edmonds-Karp implementation for finding maximum flow."""
    parent = [-1] * len(capacity)
    max_flow = 0

    while bfs(capacity, source, sink, parent):
        path_flow = float('Inf')
        s = sink

        while s != source:
            path_flow = min(path_flow, capacity[parent[s]][s])
            s = parent[s]

        max_flow += path_flow
        v = sink
        while v != source:
            u = parent[v]
            capacity[u][v] -= path_flow
            capacity[v][u] += path_flow
            v = parent[v]

    return max_flow

def find_min_cut(adj_matrix, n,total):
    """Finds a min-cut in the graph by running max-flow algorithm."""
    # Example source and sink. In a more complex case, you might want to explore different pairs.
    source = 0
    sink = sorted(n)[-1]
    last = sink
    
    capacity = adj_matrix.toarray().tolist()
    max_flow = edmonds_karp(capacity, source, sink)

    visited = [False] * total
    bfs(capacity, source, sink, visited)
    
    cut_edges = []
    for i in n:
        for j in n:
            if visited[i] and not visited[j] and adj_matrix[i, j] > 0:
                cut_edges.append((i, j))
    
    return cut_edges


def join_sccs(adj_matrix, sccs):
    """Join SCCs by adding edges between them."""
    # A simple strategy: join SCCs in sequence
    joined_edges = []
    for i in range(len(sccs) - 1):
        u = sccs[i][0]
        v = sccs[i + 1][0]
        adj_matrix[u, v] = 1  # Add a small-weight edge to connect the components
        joined_edges.append((u, v))
    
    return joined_edges

def perform_min_cut_or_join(sccs, adj_matrix, min_scc_range, max_scc_range,n):
    """Modify the graph based on SCC count."""
    length = len(sccs)
    cntter=0
    if length < min_scc_range:
        while length < min_scc_range:
            sccs1 = find_sccs(adj_matrix, n)

            for component in sccs1:
                if cntter == len(sccs1):
                    return 
                print(f"SCCs ({length}) below the min range ({min_scc_range}). Performing min-cut...")

                bridges = find_bridges(adj_matrix,component, n)
            
                if not bridges:
                    cntter+=1
                    print("No bridges left to cut.")
                    continue
                
                # Remove bridge edges
                for u, v in bridges:
                    adj_matrix[u, v] = 0
                    adj_matrix[v, u] = 0
                    # print(f"Cutting bridge between {u+1} and {v+1}")
                
                # Recalculate SCCs after cutting
                sccs = find_sccs(adj_matrix, n)
                length = len(sccs)
                if length >= min_scc_range:
                    break

    
            
    elif length > max_scc_range:

        cnt=0
        while length > max_scc_range:
            print(f"SCCs ({length}) exceed the max range ({max_scc_range}). joining...")
            u = sccs[cnt][0]
            v = sccs[cnt+1][0]
            adj_matrix[u, v] = 1
            adj_matrix[v,u] = 1
            length-=1
            cnt+=2
            if(cnt>len(sccs)):
                sccs= find_sccs(adj_matrix, n)
                length = len(sccs)

            cnt%=len(sccs)
            if cnt==len(sccs)-1:
                cnt=0
           
    else:
        print(f"SCCs ({len(sccs)}) are within the range [{min_scc_range}, {max_scc_range}]. No changes made.")
    ans = find_sccs(adj_matrix, n)

    if(ans < min_scc_range or ans> max_scc_range):
        print("Sorry range too absurd, cant cut graph more")
        
    print(f"final value: {len(ans)}")


def save_graph(filename, adj_matrix):
    """Save the graph back to the file in the same format as it was read."""
    n = adj_matrix.shape[0]
    with open(filename, 'w') as file:
        print(filename,"writing in")
        # Write the number of nodes and edges
        non_zero_elements = adj_matrix.nonzero()
        num_edges = len(non_zero_elements[0])
        file.write(f"{n} {num_edges}\n")
        
        # Write the edges with weights
        for i, j in zip(non_zero_elements[0], non_zero_elements[1]):
            # Write in the same format as the input (node1, node2, weight-1)
            weight = int(adj_matrix[i, j]) - 1
            file.write(f"{i+1} {j+1} {weight}\n")


def main(args):
    filename1 = args[0]
    filename2 = args[1]
    min_scc_range = int(args[2])
    max_scc_range = int(args[3])

    adj_matrix, n = read_graph(filename1)
    sccs = find_sccs(adj_matrix, n)
    print(len(sccs))
    perform_min_cut_or_join(sccs, adj_matrix, min_scc_range, max_scc_range,n)
    save_graph(filename1, adj_matrix)


    
    


if __name__ == "__main__":
    main(sys.argv[1:])
