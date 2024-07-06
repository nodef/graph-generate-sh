import numpy as np
import scipy.sparse as sp
import os
import matplotlib.pyplot as plt
from scipy.sparse.linalg import svds
import sys


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
        adj_matrix[node1-1, node2-1] = wt  
    
    return adj_matrix,n

def perform_svd(adj_matrix, k):
    u, s, vt = svds(adj_matrix, k=k)
    return s, u, vt

def plot_singular_values(singular_values, adj_matrix, mainname):
    # norm = np.linalg.norm(adj_matrix)
    # min_svd = np.min(singular_values)
    # cond = np.linalg.cond(adj_matrix.toarray())
    # rank = np.linalg.matrix_rank(adj_matrix.toarray())
    # null_space_dim = adj_matrix.shape[0] - rank
    # full_rank = rank == min(adj_matrix.shape)

    plt.figure(figsize=(10, 6))
    plt.scatter(range(1, len(singular_values) + 1), sorted(singular_values, reverse=True))
    plt.xlabel('Index')
    plt.ylabel('Singular Value')

    dirname = mainname.split('/')
    root = ""
    # print(mainname)
    # print(dirname)
    for i in range(len(dirname) - 1):
        root += dirname[i]
        root += "/"
    root += "svd/"
    create_folder_if_not_exists(root)
    root += dirname[-1]

    # plt.title(f'Singular Values of the Adjacency Matrix\n'
    #           f'Norm: {norm:.2e}, Min SVD: {min_svd:.2e}, Cond: {cond:.2e}, '
    #           f'Rank: {rank}, Null Space Dim: {null_space_dim}, Full Rank: {full_rank}')
    plt.grid(True)
    # print(root)
    plt.savefig(root)

def main(args):
    filename = args
    adj_matrix, n = read_graph(filename[0])
    singular_values, u, vt = perform_svd(adj_matrix, 10)
    plot_singular_values(singular_values, adj_matrix, args[1])

if __name__ == "__main__":
    main(sys.argv[1:])
