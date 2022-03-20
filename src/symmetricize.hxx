#pragma once
#include "copy.hxx"




// SYMMETRICIZE
// ------------

template <class H, class G>
void symmetricizeTo(H& a, const G& x) {
  for (int u : x.vertices())
    a.addVertex(u); // , x.vertexData(u));
  for (int u : x.vertices()) {
    for (int v : x.edges(u)) {
      a.addEdge(u, v); // , a.edgeData(u, v));
      a.addEdge(v, u); // , a.edgeData(u, v));
    }
  }
  a.correct();
}

template <class G>
auto symmetricize(const G& x) {
  G a; symmetricizeTo(a, x);
  return a;
}
