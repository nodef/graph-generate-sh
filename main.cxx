#include <tuple>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "src/main.hxx"
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




// RUN-MTX
// -------

void runMtxRewrite(const Options& o) {
  if (o.output.empty()) printf("Rewriting %s ...\n", o.file.c_str());
  else printf("Rewriting %s to %s ...\n", o.file.c_str(), o.output.c_str());
  if (o.output.empty()) rewriteMtx(cout, o.file.c_str());
  else rewriteMtx(o.output.c_str(), o.file.c_str());
}

void runMtxSamples(const Options& o) {
  printf("Loading graph %s ...\n", o.file.c_str());
  auto x  = readMtx(o.file.c_str()); println(x);
  selfLoopTo(x, [&](int u) { return isDeadEnd(x, u); });
  print(x); printf(" (selfLoopDeadEnds)\n");
  auto xt = transposeWithDegree(x);
  print(xt); printf(" (transposeWithDegree)\n");
  for (int c=0; c<o.count; c++) {
    if (o.output.empty()) printf("# DELTA %d of %d\n", c, o.count);
    else printf("Writing delta %d to %s ...\n", c, o.output.c_str());
    GraphDelta d = createMixedGraphDelta(x, o.samples/2, o.samples/2);
    if (o.output.empty()) printf("%s", toString(d).c_str());
    else toFiles(o.output.c_str(), c, d);
  }
}

void runMtx(const Options& o) {
  if (o.samples == 0) runMtxRewrite(o);
  else runMtxSamples(o);
}




// MAIN
// ----

int main(int argc, char **argv) {
  typedef FileFormat F;
  Options o = readOptions(argc, argv);
  if (o.help)           { printf("%s\n\n", helpMessage()); return 0; }
  if (!o.error.empty()) { printf("error: %s\n\n%s\n\n", o.error.c_str(), helpMessage()); return 1; }
  switch (o.format) {
    default: break;
    case F::FIXED_MTX:    runMtx(o);  break;
  }
  printf("\n");
  return 0;
}
