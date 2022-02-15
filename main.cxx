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




// RUN-MTX
// -------

void runMtx(const Options& o) {
  printf("Loading graph %s ...\n", o.file.c_str());
  auto x  = readMtx(o.file.c_str()); println(x);
  selfLoopTo(x, [&](int u) { return isDeadEnd(x, u); });
  print(x); printf(" (selfLoopDeadEnds)\n");
  auto xt = transposeWithDegree(x);
  print(xt); printf(" (transposeWithDegree)\n");
  GraphDelta d = createMixedGraphDelta(x, o.samples/2, o.samples/2);
  printf("%s", toString(d).c_str());
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
