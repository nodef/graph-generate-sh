#pragma once
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include "_main.hxx"
#include "DiGraph.hxx"

using std::string;
using std::istream;
using std::stringstream;
using std::ofstream;
using std::getline;




// PROCESS-EDGELIST
// ----------------

template <class FE, class FC>
void processEdgelist(istream& s, FE fe, FC fc) {
  string ln;
  while (getline(s, ln)) {
    int u, v;
    stringstream ls(ln);
    if (!(ls >> u >> v)) break;
    fe(u, v);
  }
  fc();
}

template <class FE, class FC>
void processEdgelist(const char *pth, FE fe, FC fc) {
  string buf = readFile(pth);
  stringstream s(buf);
  processEdgelist(s, fe, fc);
}




// READ-EDGELIST
// -------------

template <class G>
void readEdgelist(G& a, istream& s) {
  auto fe = [&](int u, int v) { a.addEdge(u, v); };
  auto fc = [&]() { a.correct(); };
  processEdgelist(s, fe, fc);
}

auto readEdgelist(istream& s) {
  DiGraph<> a; readEdgelist(a, s);
  return a;
}


template <class G>
void readEdgelist(G& a, const char *pth) {
  string buf = readFile(pth);
  stringstream s(buf);
  return readEdgelist(a, s);
}

auto readEdgelist(const char *pth) {
  DiGraph<> a; readEdgelist(a, pth);
  return a;
}




// WRITE-EDGELIST
// --------------

template <class G>
void writeEdgelist(ostream& a, const G& x) {
  for (int u : x.vertices()) {
    for (int v : x.edges(u))
      a << u << " " << v << "\n";
  }
}

template <class G>
void writeEdgelist(string pth, const G& x) {
  string s0; stringstream s(s0);
  writeEdgelist(s, x);
  ofstream f(pth);
  f << s.rdbuf();
  f.close();
}
