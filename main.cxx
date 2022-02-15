#include <tuple>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "src/main.hxx"
#include "options.hxx"

using namespace std;




// SHOW
// ----

auto groupDetails(const vector2d<int>& gs) {
  int gmin = gs.empty()? 0 : gs[0].size();
  int gmax = 0, GS = gs.size(), GE = 0;
  for (const auto& g : gs) {
    int G = g.size();
    if (G<gmin) gmin = G;
    if (G>gmax) gmax = G;
    GE += G;
  }
  return make_tuple(GS, gmin, gmax, GS>0? GE/float(GS) : 0);
}


template <class G>
void showBasics(const char *pre, const G& x) {
  printf("- %sorder:   %d (vertices)\n",        pre, x.order());
  printf("- %ssize:    %d (edges)\n",           pre, x.size());
  printf("- %sspan:    %d (vertex-id limit)\n", pre, x.span());
  printf("- %sdensity: %e (fully connectedness fraction)\n", pre, density(x));
}

template <class G, class H>
void showDegree(const char *pre, const G& x, const H& xt) {
  printf("- %sdegree-min:   %d (minimum out-degree)\n", pre, minDegree(x));
  printf("- %sdegree-max:   %d (maximum out-degree)\n", pre, maxDegree(x));
  printf("- %sdegree-avg:   %f (average out-degree)\n", pre, avgDegree(x));
  printf("- %sindegree-min: %d (minimum in-degree)\n",  pre, minDegree(xt));
  printf("- %sindegree-max: %d (maximum in-degree)\n",  pre, maxDegree(xt));
  printf("- %sindegree-avg: %f (average in-degree)\n",  pre, avgDegree(xt));
}

template <class G, class H>
void showSpecial(const char *pre, const G& x, const H& xt) {
  printf("- %sdead-ends:      %d (vertices with no outgoing edges)\n", pre, deadEndCount(x));
  printf("- %sdead-ends-rec:  %d (vertices pointing to dead ends)\n",  pre, recursiveDeadEndCount(x));
  printf("- %sbeginnings:     %d (vertices with no incoming edges)\n", pre, deadEndCount(xt));
  printf("- %sbeginnings-rec: %d (vertices pointed by beginnings)\n",  pre, recursiveDeadEndCount(xt));
  printf("- %sself-loops:     %d (vertices pointing to self)\n",       pre, selfLoopCount(x));
}

template <class G, class H>
void showLevels(const char *pre, const G& x, const H& xt) {
  auto [lnum, lmin, lmax, lavg] = groupDetails(levelwiseGroups(x, xt));
  printf("- %slevels:    %d (vertices in a level are independent)\n", pre, lnum);
  printf("- %slevel-min: %d (minimum vertices in a level)\n",         pre, lmin);
  printf("- %slevel-max: %d (maximum vertices in a level)\n",         pre, lmax);
  printf("- %slevel-avg: %f (average vertices in a level)\n",         pre, lavg);
}

template <class G, class H>
void showComponents(const char *pre, const G& x, const H& xt) {
  auto [cnum, cmin, cmax, cavg] = groupDetails(components(x, xt));
  printf("- %scomponents:    %d (strongly connected)\n",     pre, cnum);
  printf("- %scomponent-min: %d (minimum component size)\n", pre, cmin);
  printf("- %scomponent-max: %d (maximum component size)\n", pre, cmax);
  printf("- %scomponent-avg: %f (average component size)\n", pre, cavg);
}

template <class G, class H>
void showInIdenticals(const char *pre, const G& x, const H& xt) {
  auto [inum, imin, imax, iavg] = groupDetails(edgeIdenticals(xt));
  printf("- %sinidenticals:    %d (vertices with identical in-edges)\n", pre, inum);
  printf("- %sinidentical-min: %d (minimum in-identical size)\n", pre, imin);
  printf("- %sinidentical-max: %d (maximum in-identical size)\n", pre, imax);
  printf("- %sinidentical-avg: %f (average in-identical size)\n", pre, iavg);
}

template <class G, class H>
void showOutIdenticals(const char *pre, const G& x, const H& xt) {
  auto [onum, omin, omax, oavg] = groupDetails(edgeIdenticals(x));
  printf("- %soutidenticals:    %d (vertices with identical out-edges)\n", pre, onum);
  printf("- %soutidentical-min: %d (minimum out-identical size)\n", pre, omin);
  printf("- %soutidentical-max: %d (maximum out-identical size)\n", pre, omax);
  printf("- %soutidentical-avg: %f (average out-identical size)\n", pre, oavg);
}

template <class G, class H>
void showChains(const char *pre, const G& x, const H& xt) {
  auto [cnum, cmin, cmax, cavg] = groupDetails(chains(xt, x));
  printf("- %schains:    %d (vertices connected in a line)\n", pre, cnum);
  printf("- %schain-min: %d (minimum chain size)\n", pre, cmin);
  printf("- %schain-max: %d (maximum chain size)\n", pre, cmax);
  printf("- %schain-avg: %f (average chain size)\n", pre, cavg);
}

template <class G, class H>
void showBlockgraph(const char *pre, const G& x, const H& xt) {
  auto cs = components(x, xt);
  auto b  = blockgraph(x, cs);
  auto bt = transposeWithDegree(b);
  string sf(pre); sf += "blockgraph-";
  const char *prf = sf.c_str();
  showBasics (prf, b);
  showDegree (prf, b, bt);
  showSpecial(prf, b, bt);
  showLevels (prf, b, bt);
}


template <class G, class H>
void showAll(const char *pre, const G& x, const H& xt, const Options& o) {
  showBasics (pre, x);
  showDegree (pre, x, xt);
  showSpecial(pre, x, xt);
  // showDepth  (pre, x, xt);  // slow!
  if (o.components) showComponents   (pre, x, xt);
  if (o.identicals) showInIdenticals (pre, x, xt);
  if (o.identicals) showOutIdenticals(pre, x, xt);
  if (o.chains)     showChains       (pre, x, xt);
  if (o.blockgraph) showBlockgraph   (pre, x, xt);
}




// DO TRANSFORM
// ------------

template <class G>
void doTransformTo(G& a, GraphTransform o) {
  typedef GraphTransform T;
  switch (o) {
    default: break;
    case T::IDENTITY: break;
    case T::LOOP_DEADENDS:
      selfLoopTo(a, [&](int u) { return isDeadEnd(a, u); });
      print(a); printf(" (loopDeadEnds)\n"); break;
    case T::LOOP_VERTICES:
      selfLoopTo(a, [&](int u) { return true; });
      print(a); printf(" (loopVertices)\n"); break;
  }
}

template <class G>
auto doTransform(const G& x, GraphTransform o) {
  auto a = duplicate(x); doTransformTo(a, o);
  return a;
}




// RUN-MTX
// -------

void runMtx(const Options& o) {
  printf("Loading graph %s ...\n", o.file.c_str());
  auto x  = readMtxOutDiGraph(o.file.c_str()); println(x);
  doTransformTo(x, o.transform);
  auto xt = transposeWithDegree(x);
  print(xt); printf(" (transposeWithDegree)\n");
  showAll("", x, xt, o);
}




// RUN-SNAP
// --------

void runSnap(const Options& o) {
  printf("Using graph %s ...\n", o.file.c_str());
  string data = readFile(o.file.c_str());
  int M = countLines(data), steps = o.samples, jump = M/steps;
  printf("- temporal-edges: %d (including repeated ones)\n", M);
  printf("---\n");
  DiGraph<> xo;
  stringstream s(data);
  for (int i=0;; ++i) {
    if (!readSnapTemporalTo(xo, s, jump)) break; println(xo);
    auto x  = doTransform(xo, o.transform);
    auto xt = transposeWithDegree(x);
    print(xt); printf(" (transposeWithDegree)\n");
    printf("- sample:  %d (of %d)\n", i+1, o.samples);
    showAll("", x, xt, o);
    printf("---\n");
  }
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
    case F::TEMPORAL_TXT: runSnap(o); break;
  }
  printf("\n");
  return 0;
}
