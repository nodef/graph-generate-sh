import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.sparse.linalg import svds

def read_graph(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    n, m = map(int, lines[0].split())
    
    adj_matrix = sp.lil_matrix((n, n), dtype=np.float64)
    
    for line in lines[1:]:
        node1, node2, wt = map(int, line.split())
        adj_matrix[node1-1, node2-1] = wt  
    
    return adj_matrix,n

def perform_svd(adj_matrix, k):
    u, s, vt = svds(adj_matrix, k=k)
    return s

def plot_singular_values(singular_values):
    plt.figure(figsize=(10, 6))
    plt.scatter(range(1, len(singular_values) + 1), sorted(singular_values, reverse=True))
    plt.xlabel('Index')
    plt.ylabel('Singular Value')
    plt.title('Singular Values of the Adjacency Matrix')
    plt.grid(True)
    # plt.show()
    plt.savefig()

def main():
    filename = 'output/custom_1'
    adj_matrix,n = read_graph(filename)
    singular_values = perform_svd(adj_matrix,100)
    plot_singular_values(singular_values)

if __name__ == "__main__":
    main()
