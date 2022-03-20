#pragma once
#include "copy.hxx"




// SYMMETRICIZE
// ------------

template <class G>
void symmetricizeTo(G& a) {
  for (int u : a.vertices()) {
    for (int v : a.edges(u))
      a.addEdge(v, u, a.edgeData(u, v));
  }
  a.correct();
}

template <class G>
auto symmetricize(const G& x) {
  auto a = copy(x); symmetricizeTo(a);
  return a;
}
