#include <tuple>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "_main.hxx"
#include "update.hxx"

using std::tuple;
using std::vector;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::discrete_distribution;
using std::make_tuple;
using std::sort;
using std::unique;
using std::remove_if;
using std::exp;
using std::log;




#pragma region METHODS
#pragma region REMOVE RANDOM EDGE
/**
 * Randomly decide upon a edge that may be removed from a graph, from a given source vertex.
 * @param rnd random number generator (updated)
 * @param x input graph
 * @param u source vertex
 * @param fe callback function to remove edge (u, v, w)
 * @returns true if edge was removed
 */
template <class R, class G, class K, class FE>
inline bool removeRandomEdgeFrom(R& rnd, const G& x, K u, FE fe) {
  uniform_real_distribution<> dis(0.0, 1.0);
  if (x.degree(u) == 0) return false;
  K vi = K(dis(rnd) * x.degree(u)), i = 0;
  bool removed = false, skip = false;
  x.forEachEdge(u, [&](auto v, auto w) {
    if (skip) return;
    if (i++ == vi) { removed = fe(u, v, w); skip = true; }
  });
  return removed;
}


/**
 * Randomly decide upon a edge that may be removed from a graph.
 * @param rnd random number generator (updated)
 * @param x input graph
 * @param i begin vertex
 * @param n number of vertices / range
 * @param fe callback function to remove edge (u, v, w)
 * @returns true if edge was removed
 */
template <class R, class G, class FE>
inline bool removeRandomEdge(R& rnd, const G& x, size_t i, size_t n, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  return removeRandomEdgeFrom(rnd, x, u, fe);
}
#pragma endregion




#pragma region ADD RANDOM EDGE
/**
 * Randomly decide upon a edge that may be added to a graph.
 * @param rnd random number generator (updated)
 * @param x input graph
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param w edge weight
 * @param fe callback function to add edge (u, v, w)
 * @returns true if edge was added
 */
template <class R, class G, class V, class FE>
inline bool addRandomEdge(R& rnd, const G& x, size_t i, size_t n, V w, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  K v = K(i + n*dis(rnd));
  return fe(u, v, w);
}
#pragma endregion




#pragma region GENERATE BATCH
/**
 * Generate a batch of random edge deletions.
 * @param rnd random number generator (updated)
 * @param x original graph
 * @param batchSize number of edges to remove
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param undirected remove undirected edges?
 * @returns deleted edges {u, v}
 */
template <class R, class G>
inline auto generateEdgeDeletions(R& rnd, const G& x, size_t batchSize, size_t i, size_t n, bool undirected) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  int retries = 5;
  vector<tuple<K, K, V>> deletions;
  auto fe = [&](auto u, auto v, auto w) {
    deletions.push_back(make_tuple(u, v, w));
    if (undirected) deletions.push_back(make_tuple(v, u, w));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(rnd, x, i, n, fe); }, retries);
  return deletions;
}


/**
 * Generate a batch of random edge insertions.
 * @param rnd random number generator (updated)
 * @param x original graph
 * @param batchSize number of edges to add
 * @param i begin vertex
 * @param n number of vertices / range
 * @param undirected add undirected edges?
 * @param w edge weight
 * @returns inserted edges {u, v, w}
 */
template <class R, class G, class V>
inline auto generateEdgeInsertions(R& rnd, const G& x, size_t batchSize, size_t i, size_t n, bool undirected, V w) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K, V>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    insertions.push_back(make_tuple(u, v, w));
    if (undirected) insertions.push_back(make_tuple(v, u, w));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(rnd, x, i, n, w, fe); }, retries);
  return insertions;
}
#pragma endregion




#pragma region TIDY BATCH
/**
 * Filter out edges in batch update by existence.
 * @param edges edges in batch update (updated)
 * @param x original graph
 * @param exists true to keep existing edges, false to keep non-existing edges
 */
template <class G, class K, class V>
inline void filterEdgesByExistenceU(vector<tuple<K, K, V>>& edges, const G& x, bool exists) {
  auto ft = [&](const auto& e) {
    auto [u, v, w] = e;
    return x.hasEdge(u, v) != exists;
  };
  auto it = remove_if(edges.begin(), edges.end(), ft);
  edges.erase(it, edges.end());
}


/**
 * Sort edges in batch update by source/destination vertex.
 * @param edges edges in batch update (updated)
 */
template <class K, class V>
inline void sortEdgesByIdU(vector<tuple<K, K, V>>& edges) {
  auto fl = [](const auto& a, const auto& b) {
    auto [u1, v1, w1] = a;
    auto [u2, v2, w2] = b;
    return u1 < u2 || (u1 == u2 && v1 < v2);
  };
  sort(edges.begin(), edges.end(), fl);
}


/**
 * Keep only unique edges in batch update.
 * @param edges edges in batch update (updated)
 */
template <class K, class V>
inline void uniqueEdgesU(vector<tuple<K, K, V>>& edges) {
  auto fe = [](const auto& a, const auto& b) {
    auto [u1, v1, w1] = a;
    auto [u2, v2, w2] = b;
    return u1 == u2 && v1 == v2;
  };
  auto it = unique(edges.begin(), edges.end(), fe);
  edges.erase(it, edges.end());
}


/**
 * Filter out edges in batch update by existence, sort by source/destination vertex, and keep only unique edges.
 * @param deletions edge deletions in batch update (updated)
 * @param insertions edge insertions in batch update (updated)
 * @param x original graph
 */
template <class G, class K, class V>
inline void tidyBatchUpdateU(vector<tuple<K, K, V>>& deletions, vector<tuple<K, K, V>>& insertions, const G& x) {
  filterEdgesByExistenceU(deletions,  x, true);
  filterEdgesByExistenceU(insertions, x, false);
  sortEdgesByIdU(deletions);
  sortEdgesByIdU(insertions);
  uniqueEdgesU(deletions);
  uniqueEdgesU(insertions);
}
#pragma endregion




#pragma region APPLY
/**
 * Apply a batch update to a graph.
 * @param a input graph (updated)
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 */
template <class G, class K, class V>
inline void applyBatchUpdateU(G& a, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions) {
  for (auto [u, v, w] : deletions)
    a.removeEdge(u, v);
  updateU(a);
  for (auto [u, v, w] : insertions)
    a.addEdge(u, v, w);
  updateU(a);
}


#ifdef OPENMP
/**
 * Apply a batch update to a graph.
 * @param a input graph (updated)
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 */
template <class G, class K, class V>
inline void applyBatchUpdateOmpU(G& a, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions) {
  for (auto [u, v, w] : deletions)
    a.removeEdge(u, v);
  updateOmpU(a);
  for (auto [u, v, w] : insertions)
    a.addEdge(u, v, w);
  updateOmpU(a);
}
#endif
#pragma endregion
#pragma endregion


#pragma region UNIFORM

/**
 * Uniformly update a graph.
 * @param rng random number generator 
 * @param graph graph to update
 * @param batchSize number of edges to update
 * @param edgeInsertions fraction of edges to insert
 * @param edgeDeletions fraction of edges to delete
 * @param insertions edge insertions in batch update 
 * @param deletions edge deletions in batch update 
 * @param allowDuplicateEdges allow duplicate edges in batch?
*/
template <class R, class G>
void uniformUpdate(R& rng, G& graph, size_t batchSize, double edgeInsertions, double edgeDeletions, std::vector<std::tuple<typename G::key_type, typename G::key_type, typename G::edge_value_type>>& insertions, std::vector<std::tuple<typename G::key_type, typename G::key_type, typename G::edge_value_type>>& deletions, bool allowDuplicateEdges) {
    using K = typename G::key_type;
    using V = typename G::edge_value_type;
    size_t numDeletions = static_cast<size_t>(batchSize * edgeDeletions);
    deletions = generateEdgeDeletions(rng, graph, numDeletions, 1, graph.order(), false);
    size_t numInsertions = static_cast<size_t>(batchSize * edgeInsertions);
    insertions = generateEdgeInsertions(rng, graph, numInsertions, 1, graph.order(), false, V());
    if (!allowDuplicateEdges) tidyBatchUpdateU(deletions, insertions, graph);
}

#pragma endregion


#pragma region PREFERENTIAL

/** 
 * Preferential attachment/detachment to update a graph.
 * @param rng random number generator 
 * @param graph graph to update
 * @param batchSize number of edges to update
 * @param edgeInsertions fraction of edges to insert
 * @param edgeDeletions fraction of edges to delete
 * @param insertions edge insertions in batch update 
 * @param deletions edge deletions in batch update 
 * @param allowDuplicateEdges allow duplicate edges in batch?

*/
template <class R, class G>
void preferentialUpdate(R& rng, G& graph, size_t batchSize, double edgeInsertions, double edgeDeletions, vector<tuple<typename G::key_type, typename G::key_type, typename G::edge_value_type>>& insertions, vector<tuple<typename G::key_type, typename G::key_type, typename G::edge_value_type>>& deletions, bool allowDuplicateEdges) {
    using K = typename G::key_type;
    using V = typename G::edge_value_type;
    size_t numDeletions = static_cast<size_t>(batchSize * edgeDeletions);
    size_t numInsertions = static_cast<size_t>(batchSize * edgeInsertions);
    auto getInDegree = [&](K u) {
        size_t degree = 0;
        graph.forEachVertexKey([&](K v) {if (graph.hasEdge(v, u)) degree++;});
        return degree;
    };
    double beta = 1.0;
    double lambda = 0.0;
    vector<double> inDegreeWeights;
    double totalInDegree = 0;
    graph.forEachVertexKey([&](K u) {
        double degree = getInDegree(u);
        inDegreeWeights.push_back(exp(beta * log(1 + degree)) - lambda);
        totalInDegree += degree;
    });
    uniform_int_distribution<K> sourceDistribution(1, graph.order());
    discrete_distribution<K> targetDistribution(inDegreeWeights.begin(), inDegreeWeights.end());
    for (size_t i = 0; i < numDeletions; i++) {
        K u = sourceDistribution(rng);
        K v = targetDistribution(rng)+1;
        if (graph.hasEdge(u, v)) {
            V weight = graph.edgeValue(u, v);
            deletions.emplace_back(u, v, weight);
        }
    }
    for (size_t i = 0; i < numInsertions; i++) {
        K u = sourceDistribution(rng);
        K v = targetDistribution(rng)+1;
        if (!graph.hasEdge(u, v)) insertions.emplace_back(u, v, V());
    }
  if(!allowDuplicateEdges) tidyBatchUpdateU(deletions, insertions, graph);
}
#pragma endregion