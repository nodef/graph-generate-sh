#pragma once
#include <algorithm>
#include "copy.hxx"

using std::min;
using std::max;




// UNSYMMETRICIZE
// --------------

template <class H, class G>
void unsymmetricizeTo(H& a, const G& x) {
  for (int u : x.vertices())
    a.addVertex(u);
  for (int u : x.vertices()) {
    for (int v : x.edges(u)) {
      if (u>v) {
        a.removeEdge(u, v);
      }
      else {
        a.addEdge(v, u);
        a.removeEdge(u, v);
      }
    }
  }
  a.correct();
}

template <class G>
auto unsymmetricize(const G& x) {
  G a;
  for (int u : x.vertices())
    a.addVertex(a);
  for (int u : x.vertices()) {
    for (int v : x.edges(u)) {
      int _u = min(u, v);
      int _v = max(u, v);
      a.addEdge(_u, _v);
    }
  }
  return a;
}




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
