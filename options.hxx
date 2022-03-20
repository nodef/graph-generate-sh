#pragma once
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::stoi;




// COMMAND
// -------

enum class Command {
  UNKNOWN,
  DELTA,
  REWRITE,
};

auto parseCommand(const string& x) {
  typedef Command C;
  if (x=="delta"   || x=="batch")  return C::DELTA;
  if (x=="rewrite" || x=="recast") return C::REWRITE;
  return C::UNKNOWN;
}




// FILE FORMAT
// -----------

enum class FileFormat {
  UNKNOWN,
  FIXED_EDGES,
  FIXED_MTX,
  TEMPORAL_TXT
};

auto parseFileFormat(const string& x) {
  typedef FileFormat F;
  if (x=="fixed-edgelist" || x=="edgelist"      || x=="edges" || x==".edges") return F::FIXED_EDGES;
  if (x=="fixed-mtx"      || x=="matrix-market" || x=="mtx"   || x==".mtx")   return F::FIXED_MTX;
  if (x=="temporal-txt"   || x=="temporal"      || x=="txt"   || x==".txt")   return F::TEMPORAL_TXT;
  return F::UNKNOWN;
}

auto parseFileFormats(const string& x) {
  typedef FileFormat F; vector<F> a;
  string y = removeAll(x, ' ');
  for (size_t i=0; i<y.length();) {
    size_t p = y.find_first_of(',', i);
    if (p==string::npos) p = y.length();
    string f = y.substr(i, p-i);
    a.push_back(parseFileFormat(f));
    i = p+1;
  }
  return a;
}




// GRAPH TRANSFORM
// ---------------

enum class GraphTransform {
  UNKNOWN,
  SYMMETRICIZE,
  LOOP_DEADENDS,
  LOOP_VERTICES,
  CLEAR_WEIGHTS,
  SET_WEIGHTS,
};

auto parseGraphTransform(const string& x) {
  typedef GraphTransform T;
  if (x=="symmetricize"  || x=="symmetric"    || x=="sym") return T::SYMMETRICIZE;
  if (x=="loop-deadends" || x=="loop")     return T::LOOP_DEADENDS;
  if (x=="loop-vertices" || x=="loop-all") return T::LOOP_VERTICES;
  if (x=="clear-weights" || x=="zero-weights" || x=="no-weights"     || x=="weights=0") return T::CLEAR_WEIGHTS;
  if (x=="set-weights"   || x=="unit-weights" || x=="common-weights" || x=="weights=1") return T::SET_WEIGHTS;
  return T::UNKNOWN;
}

auto parseGraphTransforms(const string& x) {
  typedef GraphTransform T; vector<T> a;
  string y = removeAll(x, ' ');
  for (size_t i=0; i<y.length();) {
    size_t p = y.find_first_of(',', i);
    if (p==string::npos) p = y.length();
    string t = y.substr(i, p-i);
    a.push_back(parseGraphTransform(t));
    i = p+1;
  }
  return a;
}




// OPTIONS
// -------

struct Options {
  private:
  typedef Command        C;
  typedef FileFormat     F;
  typedef GraphTransform T;
  public:
  bool   help    = false;
  string error   = "";
  string input   = "";
  string output  = "";
  string commandStr    = "";
  string formatsStr    = "";
  string transformsStr = "";
  string samplesStr    = "";
  string countStr      = "";
  vector<T> transforms {};
  vector<F> formats {};
  C command   = C::UNKNOWN;
  int samples = 0;
  int count   = 1;
};




string pathExtname(const string& path) {
  auto idx = path.rfind('.');
  return idx==string::npos? "" : path.substr(idx);
}


Options readOptions(int argc, char **argv) {
  typedef Command        C;
  typedef FileFormat     F;
  typedef GraphTransform T;
  Options a;
  for (int i=1; i<argc; ++i) {
    string k = argv[i];
    if (k=="--help") a.help = true;
    else if (k=="-f" || k=="--format"    || k=="--formats")    a.formatsStr    = argv[++i];
    else if (k=="-t" || k=="--transform" || k=="--transforms") a.transformsStr = argv[++i];
    else if (k=="-s" || k=="--samples") a.samplesStr = argv[++i];
    else if (k=="-c" || k=="--count")   a.countStr   = argv[++i];
    else if (k=="-i" || k=="--input")   a.input  = argv[++i];
    else if (k=="-o" || k=="--output")  a.output = argv[++i];
    else if (k[0]=='-')            { a.error      = "\'"+k+"\' is not an option";          break; }
    else if (a.commandStr.empty()) { a.commandStr = argv[i]; }
    else if (!a.input.empty())     { a.error      = "\'"+k+"\' file cannot be read as well"; break; }
    else a.input = argv[i];
  }
  a.command = parseCommand(a.commandStr);
  if (a.command==C::UNKNOWN) { a.error = "no command specified"; return a; }
  if (a.input.empty()) { a.error = "no input file specified"; return a; }
  a.formats = parseFileFormats(a.formatsStr);
  for (F f : a.formats)
    if (f==F::UNKNOWN) { a.error = "\'"+a.formatsStr +"\' format is not recognized"; return a; }
  if (a.formats.size()<1) a.formats.push_back(parseFileFormat(pathExtname(a.input)));
  if (a.formats.size()<2) a.formats.push_back(parseFileFormat(pathExtname(a.output)));
  if (a.formats[0]==F::UNKNOWN) { a.error = "unknown input format"; return a; }
  if (a.formats[1]==F::UNKNOWN && !a.output.empty()) { a.error = "unknown output format"; return a; }
  a.transforms = parseGraphTransforms(a.transformsStr);
  for (T t : a.transforms)
    if (t==T::UNKNOWN) { a.error = "\'"+a.transformsStr +"\' transform is not recognized"; return a; }
  try { if (!a.samplesStr.empty()) a.samples = stoi(a.samplesStr); }
  catch (...) { a.error = "\'"+a.samplesStr+"\' samples is not an integer"; return a; }
  try { if (!a.countStr.empty()) a.count = stoi(a.countStr); }
  catch (...) { a.error = "\'"+a.countStr+"\' count is not an integer"; return a; }
  if (a.samples<0) { a.error = "\'"+a.samplesStr+  "\' samples must be positive"; return a; }
  if (a.count<0)   { a.error = "\'"+a.countStr+    "\' count must be positive"; return a; }
  return a;
}




// HELP
// ----

const char* helpMessage() {
  return "For usage details, please try the following URL:\n"
  "https://github.com/puzzlef/graph-generate";
}
