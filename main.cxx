#include <tuple>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "inc/main.hxx"
#include "options.hxx"

using namespace std;




// GRAPH-DELTA
// -----------

struct GraphDelta {
  vector<pair<int, int>> deletions;
  vector<pair<int, int>> insertions;
};


template <class G>
auto createMixedGraphDelta(const G& x, int del, int ins) {
  GraphDelta a;
  random_device dev;
  default_random_engine rnd(dev());
  for (int i=0; i<del; ++i)
    a.deletions.push_back(suggestRemoveRandomEdgeByDegree(x, rnd));
  for (int i=0; i<ins; ++i)
    a.insertions.push_back(suggestAddRandomEdgeByDegree(x, rnd, x.span()));
  return a;
}


string toString(const GraphDelta& x) {
  string s = ""; stringstream a(s);
  for (auto [u, v] : x.deletions)
    a << "- " << u << " " << v << "\n";
  for (auto [u, v] : x.insertions)
    a << "+ " << u << " " << v << "\n";
  return a.str();
}

void toFiles(const char* dpth, const char *ipth, const GraphDelta& x) {
  string ds; stringstream d(ds);
  string is; stringstream i(is);
  writePlain(d, x.deletions);
  writePlain(i, x.insertions);
  writeFile(dpth, d.str());
  writeFile(ipth, i.str());
}
void toFiles(const char* pth, int i, const GraphDelta& x) {
  string fil = pth;
  size_t e   = fil.find_last_of('.');
  string nam = e != size_t(-1)? fil.substr(0, e) : nam;
  string ext = e != size_t(-1)? fil.substr(e, fil.length() - e) : ".txt";
  string dstr; stringstream ds(dstr);
  string istr; stringstream is(istr);
  ds << nam << "-" << i << ext;
  is << nam << "+" << i << ext;
  toFiles(ds.str().c_str(), is.str().c_str(), x);
}




// RUN
// ---

template <class G>
void performTransform(G& a, int& w, string t) {
  if (t=="UNKNOWN") {}
  else if (t[0]=='-') {
    bool sym = t[1]=='-';
    string f = readFile(t.c_str()+(sym? 2:1));
    stringstream fs(f);
    auto fe = [&](int u, int v) { a.removeEdge(u, v); };
    auto fc = [&]() { a.correct(); };
    processEdges(fs, fe, fc, sym);
    print(a); printf(" (%s)\n", t.c_str());
  }
  else if (t[0]=='+') {
    bool sym = t[1]=='+';
    string f = readFile(t.c_str()+(sym? 2:1));
    stringstream fs(f);
    auto fe = [&](int u, int v) { a.addEdge(u, v); };
    auto fc = [&] { a.correct(); };
    processEdges(fs, fe, fc, sym);
    print(a); printf(" (%s)\n", t.c_str());
  }
  else if (t=="TRANSPOSE") {
    a = transpose(a);
    print(a); printf(" (transpose)\n");
  }
  else if (t=="UNSYMMETRICIZE") {
    a = unsymmetricize(a);
    print(a); printf(" (unsymmetricize)\n");
  }
  else if (t=="SYMMETRICIZE") {
    a = symmetricize(a);
    print(a); printf(" (symmetricize)\n");
  }
  else if (t=="LOOP_VERTICES") {
    selfLoopTo(a, [](int u) { return u; });
    print(a); printf(" (selfLoopVertices)\n");
  }
  else if (t=="LOOP_DEADENDS") {
    selfLoopTo(a, [&](int u) { return isDeadEnd(a, u); });
    print(a); printf(" (selfLoopDeadEnds)\n");
  }
  else if (t=="CLEAR_WEIGHTS") {
    w = 0;
    print(a); printf(" (clearWeights)\n");
  }
  else if (t=="SET_WEIGHTS") {
    w = 1;
    print(a); printf(" (setWeights)\n");
  }
}


void runDelta(const Options& o) {
  printf("Loading graph %s ...\n", o.input.c_str());
  auto x  = readMtx(o.input.c_str()); println(x); int w = -1;
  for (auto t : o.transforms)
    performTransform(x, w, t);
  // auto xt = transposeWithDegree(x);
  // print(xt); printf(" (transposeWithDegree)\n");
  for (int c=0; c<o.count; c++) {
    if (o.output.empty()) printf("# DELTA %d of %d\n", c, o.count);
    else printf("Writing delta %d to %s ...\n", c, o.output.c_str());
    GraphDelta d = createMixedGraphDelta(x, o.samples/2, o.samples/2);
    if (o.output.empty()) printf("%s", toString(d).c_str());
    else toFiles(o.output.c_str(), c, d);
  }
}


void runPlainRewrite(const Options& o) {
  if (o.output.empty()) { rewriteMtxEdgelist(cout, o.input.c_str()); return; }
  printf("Rewriting %s to %s ...\n", o.input.c_str(), o.output.c_str());
  rewriteMtxEdgelist(o.output.c_str(), o.input.c_str());
}

void runRewrite(const Options& o) {
  typedef FileFormat F;
  if (o.formats[1] == F::FIXED_EDGES && o.transforms.size() == 0) return runPlainRewrite(o);
  printf("Loading graph %s ...\n", o.input.c_str());
  auto x  = readMtx(o.input.c_str()); println(x); int w = -1;
  for (auto t : o.transforms)
    performTransform(x, w, t);
  if (o.formats[1] == F::FIXED_EDGES) {
    if (o.output.empty()) { writeEdgelist(cout, x); return; }
    printf("Writing to %s ...\n", o.output.c_str());
    writeEdgelist(o.output.c_str(), x); return;
  }
  string ws = w < 0? "" : to_string(w);
  auto   fw = [&](int u, int v) { return ws; };
  if (o.output.empty()) { writeMtx(cout, x, fw); return; }
  printf("Writing to %s ...\n", o.output.c_str());
  writeMtx(o.output, x, fw);
}




// MAIN
// ----

int main(int argc, char **argv) {
  typedef Command    C;
  typedef FileFormat F;
  Options o = readOptions(argc, argv);
  for (auto f : o.formats)
    if (f == F::TEMPORAL_TXT) o.error = "Temporal (.txt) format is not supported.";
  if (o.help)           { printf("%s\n\n", helpMessage()); return 0; }
  if (!o.error.empty()) { printf("error: %s\n\n%s\n\n", o.error.c_str(), helpMessage()); return 1; }
  if (o.command == C::DELTA) runDelta(o);
  else runRewrite(o);
  printf("\n");
  return 0;
}
