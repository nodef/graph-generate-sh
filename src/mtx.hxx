#pragma once
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "_main.hxx"
#include "DiGraph.hxx"

using std::string;
using std::istream;
using std::stringstream;
using std::ofstream;
using std::getline;
using std::max;




// PROCESS-MTX
// -----------

template <class FE, class FC>
void processEdges(istream& s, FE fe, FC fc, bool sym=false) {
  string ln;
  stringstream ls(ln);
  while (getline(s, ln)) {
    int u, v;
    ls = stringstream(ln);
    if (!(ls >> u >> v)) break;
    fe(u, v);
    if (sym) fe(v, u);
  }
  fc();
}

template <class FV, class FE, class FC>
void processMtx(istream& s, FV fv, FE fe, FC fc) {
  string ln, h0, h1, h2, h3, h4;

  // read header
  while (1) {
    getline(s, ln);
    if (ln.find('%')!=0) break;
    if (ln.find("%%")!=0) continue;
    stringstream ls(ln);
    ls >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  if (h1!="matrix" || h2!="coordinate") return;
  bool sym = h4=="symmetric" || h4=="skew-symmetric";

  // read rows, cols, size
  int r, c, sz;
  stringstream ls(ln);
  ls >> r >> c >> sz;
  int n = max(r, c);
  for (int u=1; u<=n; u++)
    fv(u);

  // read edges (from, to)
  processEdges(s, fe, fc, sym);
}

template <class FV, class FE, class FC>
void processMtx(const char *pth, FV fv, FE fe, FC fc) {
  string buf = readFile(pth);
  stringstream s(buf);
  processMtx(s, fv, fe, fc);
}




// READ-MTX
// --------

template <class G>
void readMtx(G& a, istream& s) {
  auto fv = [&](int u) { a.addVertex(u); };
  auto fe = [&](int u, int v) { a.addEdge(u, v); };
  auto fc = [&]() { a.correct(); };
  processMtx(s, fv, fe, fc);
}

auto readMtx(istream& s) {
  DiGraph<> a; readMtx(a, s);
  return a;
}


template <class G>
void readMtx(G& a, const char *pth) {
  string buf = readFile(pth);
  stringstream s(buf);
  return readMtx(a, s);
}

auto readMtx(const char *pth) {
  DiGraph<> a; readMtx(a, pth);
  return a;
}




// WRITE-MTX
// ---------

template <class G, class FW>
void writeMtx(ostream& a, const G& x, FW fw) {
  a << "%%MatrixMarket matrix coordinate real asymmetric\n";
  a << x.order() << " " << x.order() << " " << x.size() << "\n";
  for (int u : x.vertices()) {
    for (int v : x.edges(u))
      a << u << " " << v << " " << fw(u, v) << "\n";
  }
}
template <class G>
void writeMtx(ostream& a, const G& x) {
  auto fw = [&](int u, int v) { return x.edgeData(u, v); };
  writeMtx(a, x, fw);
}

template <class G, class FW>
void writeMtx(string pth, const G& x, FW fw) {
  string s0; stringstream s(s0);
  writeMtx(s, x, fw);
  ofstream f(pth);
  f << s.rdbuf();
  f.close();
}

template <class G>
void writeMtx(string pth, const G& x) {
  auto fw = [&](int u, int v) { return x.edgeData(u, v); };
  writeMtx(pth, x, fw);
}




// REWRITE-MTX
// -----------

void rewriteMtxEdgelist(ostream& a, istream& s) {
  auto fv = [&](int u) {};
  auto fe = [&](int u, int v) { a << u << ' ' << v << '\n'; };
  auto fc = [&]() {};
  processMtx(s, fv, fe, fc);
}
void rewriteMtxEdgelist(ostream& a, const char *pth) {
  string buf = readFile(pth);
  stringstream s(buf);
  rewriteMtxEdgelist(a, s);
}
void rewriteMtxEdgelist(const char *out, const char *pth) {
  string buf = readFile(pth), abuf;
  stringstream is(buf), a(abuf);
  rewriteMtxEdgelist(a, is);
  writeFile(out, a.str());
}
