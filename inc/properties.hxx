#pragma once
#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <limits>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <cmath>
#include "_main.hxx"
#include "bfs.hxx"
#include "dfs.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::pow;
using std::min;
using std::max;
using std::tuple;
using std::stack;
using std::numeric_limits;
using std::queue;
using std::max_element;
using std::set;




#pragma region METHODS
#pragma region GRAPH DATA
/**
 * Obtain the vertex keys of a graph.
 * @param x given graph
 * @returns vertex keys
 */
template <class G>
inline auto vertexKeys(const G& x) {
  using  K = typename G::key_type;
  size_t N = x.order();
  vector<K> a;
  a.reserve(N);
  x.forEachVertexKey([&](auto u) { a.push_back(u); });
  return a;
}


/**
 * Obtain the vertex value of each vertex.
 * @param a vertex value of each vertex (output)
 * @param x given graph
 */
template <class G, class V>
inline void vertexValuesW(vector<V>& a, const G& x) {
  x.forEachVertex([&](auto u, auto d) { a[u] = d; });
}


/**
 * Obtain the outgoing degree of each vertex.
 * @param a degrees of each vertex (output)
 * @param x given graph
 */
template <class G, class K>
inline void degreesW(vector<K>& a, const G& x) {
  x.forEachVertexKey([&](auto u) { a[u] = x.degree(u); });
}
#pragma endregion




#pragma region EDGE WEIGHT
/**
 * Find the total outgoing edge weight of a vertex.
 * @param x given graph
 * @param u given vertex
 * @returns total outgoing weight of a vertex
 */
template <class G, class K>
inline double edgeWeight(const G& x, K u) {
  double a = 0;
  x.forEachEdge(u, [&](auto v, auto w) { a += w; });
  return a;
}


/**
 * Find the total edge weight of a graph.
 * @param x given graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeight(const G& x) {
  double a = 0;
  x.forEachVertexKey([&](auto u) { a += edgeWeight(x, u); });
  return a;
}


#ifdef OPENMP
/**
 * Find the total edge weight of a graph.
 * @param x given graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeightOmp(const G& x) {
  using K = typename G::key_type;
  double a = 0;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a += edgeWeight(x, u);
  }
  return a;
}
#endif
#pragma endregion




#pragma region MODULARITY
/**
 * Find the modularity of a community C.
 * @param cin total weight of edges within community C
 * @param ctot total weight of edges of community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 * @see https://www.youtube.com/watch?v=0zuiLBOIcsw
 */
inline double modularityCommunity(double cin, double ctot, double M, double R=1) {
  ASSERT(cin>=0 && ctot>=0 && M>0 && R>0);
  return cin/(2*M) - R*pow(ctot/(2*M), 2);
}


/**
 * Find the modularity of a set of communities.
 * @param cin total weight of edges within each community
 * @param ctot total weight of edges of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class V>
inline double modularityCommunities(const vector<V>& cin, const vector<V>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  for (size_t i=0, I=cin.size(); i<I; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}


#ifdef OPENMP
/**
 * Find the modularity of a set of communities.
 * @param cin total weight of edges within each community
 * @param ctot total weight of edges of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class V>
inline double modularityCommunitiesOmp(const vector<V>& cin, const vector<V>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  size_t C = cin.size();
  #pragma omp parallel for schedule(static) reduction(+:a)
  for (size_t i=0; i<C; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}
#endif


/**
 * Find the modularity of a graph, based on community membership function.
 * @param x given graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityBy(const G& x, FC fc, double M, double R=1) {
  using  K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> cin(S), ctot(S);
  x.forEachVertexKey([&](auto u) {
    K c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      K d = fc(v);
      if (c==d) cin[c] += w;
      ctot[c] += w;
    });
  });
  return modularityCommunities(cin, ctot, M, R);
}


#ifdef OPENMP
/**
 * Find the modularity of a graph, based on community membership function.
 * @param x given graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityByOmp(const G& x, FC fc, double M, double R=1) {
  using  K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> vin(S), vtot(S);
  vector<double> cin(S), ctot(S);
  // Compute the internal and total weight of each vertex.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      K d = fc(v);
      if (c==d) vin[u] += w;
      vtot[u] += w;
    });
  }
  // Compute the internal and total weight of each community.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = fc(u);
    #pragma omp atomic
    cin[c]  += vin[u];
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
  return modularityCommunitiesOmp(cin, ctot, M, R);
}
#endif
#pragma endregion




#pragma region DELTA MODULARITY
/**
 * Find the change in modularity when moving a vertex from community D to C.
 * @param vcout total weight of edges from vertex v to community C
 * @param vdout total weight of edges from vertex v to community D
 * @param vtot total weight of edges from vertex v
 * @param ctot total weight of edges from community C
 * @param dtot total weight of edges from community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns delta-modularity [-0.5, 1]
 * @see https://gist.github.com/wolfram77/a3c95cd94a38a100f9b075594a823928
 */
inline double deltaModularity(double vcout, double vdout, double vtot, double ctot, double dtot, double M, double R=1) {
  ASSERT(vcout>=0 && vdout>=0 && vtot>=0 && ctot>=0 && dtot>=0 && M>0 && R>0);
  return (vcout-vdout)/M - R*vtot*(vtot+ctot-dtot)/(2*M*M);
}
#pragma endregion




#pragma region COMMUNITIES
/**
 * Obtain the size of each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns size of each community
 */
template <class G, class K>
inline vector<K> communitySize(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ++a[c];
  });
  return a;
}


#ifdef OPENMP
/**
 * Obtain the size of each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns size of each community
 */
template <class G, class K>
inline vector<K> communitySizeOmp(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a(S);
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    ++a[c];
  }
  return a;
}
#endif


/**
 * Obtain the vertices belonging to each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns vertices belonging to each community
 */
template <class G, class K>
inline vector2d<K> communityVertices(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    a[c].push_back(u);
  });
  return a;
}


#ifdef OPENMP
/**
 * Obtain the vertices belonging to each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns vertices belonging to each community
 */
template <class G, class K>
inline vector2d<K> communityVerticesOmp(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  #pragma omp parallel
  {
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K c = vcom[u];
      if (belongsOmp(c)) a[c].push_back(u);
    }
  }
  return a;
}
#endif


/**
 * Obtain the community ids of vertices in a graph.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns community ids
 */
template <class G, class K>
inline vector<K> communities(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a;
  vector<char> vis(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (vis[c]) return;
    vis[c] = 1;
    a.push_back(c);
  });
  return a;
}
#pragma endregion




#pragma region DISCONNECTED COMMUNITIES
#ifdef OPENMP
/**
 * Examine if each community in a graph is disconnected (using single flag vector, BFS).
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns whether each community is disconnected
 */
template <class G, class K>
inline vector<char> communitiesDisconnectedOmp(const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  int     T = omp_get_max_threads();
  auto coms = communitySizeOmp(x, vcom);
  vector<char> a(S), vis(S);
  vector2d<K>  us (T), vs(T);
  #pragma omp parallel
  {
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      int t = omp_get_thread_num();
      K   c = vcom[u], reached = K();
      if (coms[c]==0 || !belongsOmp(c)) continue;
      auto ft = [&](auto v, auto d) { return vcom[v]==c; };
      auto fp = [&](auto v, auto d) { ++reached; };
      us[t].clear(); vs[t].clear(); us[t].push_back(u);
      bfsVisitedForEachU(vis, us[t], vs[t], x, ft, fp);
      if (reached < coms[c]) a[c] = 1;
      coms[c] = 0;
    }
  }
  return a;
}
#endif
#pragma endregion

#pragma region DENSITY
/**
 * Find the density of a graph.
 * @param x given graph
 * @returns density of the graph
 */
template <class G>
inline double density(const G& x) {
  size_t N = x.order();
  size_t E = x.size();
  return static_cast<double>(E) / (N * N);
}
#pragma endregion

#pragma region DEGREE DISTRIBUTION
/**
 * Find the out-degree distribution of a graph.
 * @param x given graph
 * @returns out-degrees of each vertex
 */
template <class G, class K>
inline vector<K> degreeDistribution(const G& x) {
  size_t mn = SIZE_MAX, mx = 0, sum = 0;
  vector<K> a(x.order()+1, 0);
  x.forEachVertexKey([&](auto u){
    x.template forEachEdgeKey<std::function<void(K)>>(u, [&](auto v) { ++a[v]; });
  });
  return a;
}
#pragma endregion

#pragma region STRONGLY CONNECTED COMPONENTS
/**
 * Find the number of strongly connected components of a graph.
 * @param x given graph
 * @returns number of strongly connected components
 */
template <class G>
inline size_t tarjanSCC(const G& x) {
  using K = typename G::key_type;
  size_t N = x.order();
  vector<K> index(N+1, -1);
  vector<K> low(N+1, -1);
  vector<bool> onstack(N+1, false);
  stack<K> st;
  size_t foundat = 1;
  size_t count = 0;
  std::function<void(K)> strongconnect = [&](K u) {
    index[u] = low[u] = foundat++;
    st.push(u);
    onstack[u] = true;
    x.forEachEdgeKey(u, [&](K v) {
      if (index[v] == -1) {
        strongconnect(v);
        low[u] = min(low[u], low[v]);
      } else if (onstack[v]) {
        low[u] = min(low[u], index[v]);
      }
    });
    if (low[u] == index[u]) {
      count++;
      vector<K> sccTemp;
      K w;
      do {
        w = st.top();
        st.pop();
        onstack[w] = false;
        sccTemp.push_back(w);
      } while (w != u);
      // cout << "SCC: ";
      // for (auto v : sccTemp) cout << v << " ";
      // cout << std::endl;
    }
  };
  x.forEachVertexKey([&](K v) {
    if (index[v] == -1) {
      strongconnect(v);
    }
  });
  return count;
}



#pragma endregion

#pragma region KULLBACK LEIBLER DIVERGENCE
/**
* @brief  Calculate the KL distance between two distributions
* @param P The vector for distribution 1.
* @param Q The vector for distribution 2.
*/
double KLDivergence(const std::vector<double>& P, const std::vector<double>& Q) {
    size_t maxSize = std::max(P.size(), Q.size());
    double divergence = 0.0;

    for (size_t i = 0; i < maxSize; ++i) {
        double pVal = (i < P.size()) ? P[i] : 0.0;
        double qVal = (i < Q.size()) ? Q[i] : 0.0;

        if (pVal != 0) {
            if (qVal == 0) {
                throw std::invalid_argument("Q[i] must be non-zero where P[i] is non-zero.");
            }
            divergence += pVal *log(pVal*1.0 / qVal);
        }
    }
    return divergence;
}


/**
* @brief  Normalises the vector
* @param values The vector to be normalised.
* @return The normalised value vector
*/


std::vector<double> normalize(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    std::vector<double> normalized;
    for (double value : values) {
        normalized.push_back(value / sum);
    }
    return normalized;
}

/**
* @brief  Calculates the indegree distribution
* @param graph The graph object for which the distribution is calculated.
* @param distribution The distribution according to the weights in the graph
*/

template <typename G, typename K>
void calculateInDegreeDistribution(const G& graph, std::map<size_t, size_t>& distribution) {
    graph.forEachVertexKey([&](K u) {
        size_t inDegree = graph.indegree(u);
        distribution[inDegree]++;
    });
}

/**
* @brief  Converts the vector distribution in to probability distribution
* @param distribution The distribution according to the weights in the graph
* @return probabilities according to the weights
*/

std::vector<double> degreeDistributionToProbability(const std::map<size_t, size_t>& distribution) {
    std::vector<double> probabilities;
    size_t totalVertices = 0;
    for (const auto& pair : distribution) {
        totalVertices += pair.second;
    }
    for (const auto& pair : distribution) {
        probabilities.push_back(static_cast<double>(pair.second) / totalVertices);
    }
    return probabilities;
}

#pragma endregion

#pragma region DIAMETER

template <class K, class V, class E>
vector<int> getDistanceFromU(const DiGraph<K, V, E>& graph, K u) {
  vector<int> distanceFromU(graph.order()+1, numeric_limits<int>::max());
  vector<bool> visited(graph.order()+1, false);
  queue<K> BFSQueue;
  visited[u] = true;
  distanceFromU[u] = 0;
  BFSQueue.push(u);
  while (!BFSQueue.empty()) {
    K top = BFSQueue.front();
    BFSQueue.pop();
    graph.forEachEdgeKey(top, [&](K v) {
      if (!visited[v]) {
        visited[v] = true;
        distanceFromU[v] = distanceFromU[top] + 1;
        BFSQueue.push(v);
      }
    });
  }
  return distanceFromU;
}

template <class K, class V, class E>
int computeEccentricity(const DiGraph<K, V, E>& graph, K u) {
  vector<int> distanceFromU = getDistanceFromU(graph, u);
  return *max_element(distanceFromU.begin()+1, distanceFromU.end());
}

template <class K, class V, class E>
vector<set<K>> computeF(const DiGraph<K, V, E>& graph, K u) {
  vector<int> distanceFromU = getDistanceFromU(graph, u);
  int eccentricity = computeEccentricity(graph, u);
  vector<set<K>> F(eccentricity + 1);
  graph.forEachVertexKey([&](K v) {
    F[distanceFromU[v]].insert(v);
  });
  return F;
}

template <class K, class V, class E>
int computeMaxEccentricity(const DiGraph<K, V, E>& graph, const set<K>& vertices) {
  int maxEccentricity = 0;
  for (K v : vertices) {
    maxEccentricity = max(maxEccentricity, computeEccentricity(graph, v));
  }
  return maxEccentricity;
}

template <class K, class V, class E>
int iFUB(const DiGraph<K, V, E>& graph, K u, int l, int k) {
  int eccentricityU = computeEccentricity(graph, u);
  int i = eccentricityU;
  int lb = max(eccentricityU, l);
  int ub = 2 * eccentricityU;
  vector<set<K>> F = computeF(graph, u);
  while (ub - lb > k) {
    int newLowerBound = max(lb, computeMaxEccentricity(graph, F[i]));
    if (newLowerBound > 2 * (i - 1)) {
      return newLowerBound;
    } else {
      lb = newLowerBound;
      ub = 2 * (i - 1);
    }
    i--;
  }
  return lb;
}

template <class K, class V, class E>
K getVertexWithMaximumDegree(const DiGraph<K, V, E>& graph) {
  K maxDegreeVertex = 0;
  int maxDegree = 0;
  graph.forEachVertexKey([&](K v) {
    int degree = graph.degree(v);
    if (degree > maxDegree) {
      maxDegreeVertex = v;
      maxDegree = degree;
    }
  });
  return maxDegreeVertex;
}

template <class K, class V, class E>
K maxDistantVertex(const DiGraph<K, V, E>& graph, K u) {
  vector<E> distanceFromU = getDistanceFromU(graph, u);
  K maxDistantVertex = u;
  E maxDistance = distanceFromU[maxDistantVertex];
  graph.forEachVertexKey([&](K v) {
    if (distanceFromU[v] > maxDistance) {
      maxDistantVertex = v;
      maxDistance = distanceFromU[v];
    }
  });
  return maxDistantVertex;
}

template <class K, class V, class E>
K midVertex(const DiGraph<K, V, E>& graph, K u, K v) {
  vector<E> distanceFromU(graph.order()+1);
  vector<K> parentVertex(graph.order()+1);
  vector<bool> visited(graph.order()+1, false);
  queue<K> BFSQueue;
  visited[u] = true;
  distanceFromU[u] = 0;
  parentVertex[u] = -1;
  BFSQueue.push(u);
  while (!BFSQueue.empty()) {
    K top = BFSQueue.front();
    BFSQueue.pop();
    graph.forEachEdgeKey(top, [&](K neighbor) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        distanceFromU[neighbor] = distanceFromU[top] + 1;
        parentVertex[neighbor] = top;
        BFSQueue.push(neighbor);
      }
    });
    if (visited[v]) {
      break;
    }
  }
  E midDistance = distanceFromU[v] / 2;
  K midVertex = v;
  while (midDistance--) {
    midVertex = parentVertex[midVertex];
  }
  return midVertex;
}

template <class K, class V, class E>
pair<E, K> fourSweep(const DiGraph<K, V, E>& graph) {
  K r1 = getVertexWithMaximumDegree(graph);
  K a1 = maxDistantVertex(graph, r1);
  K b1 = maxDistantVertex(graph, a1);
  K r2 = midVertex(graph, a1, b1);
  K a2 = maxDistantVertex(graph, r2);
  K b2 = maxDistantVertex(graph, a2);
  K u = midVertex(graph, a2, b2);
  E lb = min(computeEccentricity(graph, a1), computeEccentricity(graph, a2));
  return make_pair(lb, u);
}

template <class K, class V, class E>
int getDiameter(const DiGraph<K, V, E>& graph) {
  auto sym_graph = symmetrize(graph);
  vector<bool> visited(sym_graph.order()+1, false);
  queue<K> BFSQueue;
  K start;
  sym_graph.forEachVertexKey([&](K v) {
    start = v;
  });
  visited[start] = true;
  BFSQueue.push(start);
  while (!BFSQueue.empty()) {
    K top = BFSQueue.front();
    BFSQueue.pop();
    sym_graph.forEachEdgeKey(top, [&](K neighbor) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        BFSQueue.push(neighbor);
      }
    });
  }
  int flg=1;
  sym_graph.forEachVertexKey([&](K v) {
    if(!visited[v]) {flg=0;}
  });
  if(flg==0) {
    return -1;
  }
  pair<E, K> p = fourSweep(sym_graph);
  return iFUB(sym_graph, p.second, p.first, 0);
}

#pragma endregion

#pragma endregion
