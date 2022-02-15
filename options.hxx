#pragma once
#include <string>

using std::string;
using std::stoi;




// FILE FORMAT
// -----------

enum class FileFormat {
  UNKNOWN,
  FIXED_MTX,
  TEMPORAL_TXT
};

auto parseFileFormat(const string& x) {
  typedef FileFormat F;
  if (x=="mtx" || x==".mtx") return F::FIXED_MTX;
  if (x=="txt" || x==".txt") return F::TEMPORAL_TXT;
  return F::UNKNOWN;
}




// GRAPH TRANSFORM
// ---------------

enum class GraphTransform {
  UNKNOWN,
  IDENTITY,
  LOOP_DEADENDS,
  LOOP_VERTICES
};

auto parseGraphTransform(const string& x) {
  typedef GraphTransform T;
  if (x=="identity" || x=="default" || x=="none" || x=="") return T::IDENTITY;
  if (x=="loop-deadends" || x=="loop")     return T::LOOP_DEADENDS;
  if (x=="loop-vertices" || x=="loop-all") return T::LOOP_VERTICES;
  return T::UNKNOWN;
}




// OPTIONS
// -------

struct Options {
  private:
  typedef FileFormat     F;
  typedef GraphTransform T;
  public:
  bool   help  = false;
  string error = "";
  string file  = "";
  string formatStr    = "";
  string transformStr = "";
  string samplesStr   = "";
  F format    = F::UNKNOWN;
  T transform = T::IDENTITY;
  int  samples    = 10;
  bool components = true;
  bool blockgraph = true;
  bool chains     = true;
  bool identicals = true;
};




string pathExtname(const string& path) {
  auto idx = path.rfind('.');
  return idx==string::npos? "" : path.substr(idx);
}


Options readOptions(int argc, char **argv) {
  typedef FileFormat     F;
  typedef GraphTransform T;
  Options a;
  for (int i=1; i<argc; ++i) {
    string k = argv[i];
    if (k=="--help") a.help = true;
    else if (k=="-f" || k=="--format")    a.formatStr    = argv[++i];
    else if (k=="-t" || k=="--transform") a.transformStr = argv[++i];
    else if (k=="-s" || k=="--samples")   a.samplesStr   = argv[++i];
    else if (k=="--components") a.components = true;
    else if (k=="--blockgraph") a.blockgraph = true;
    else if (k=="--chains")     a.chains     = true;
    else if (k=="--identicals") a.identicals = true;
    else if (k[0]=='-')       { a.error = "\'"+k+"\' is not an option";            break; }
    else if (!a.file.empty()) { a.error = "\'"+k+"\' file cannot be read as well"; break; }
    else a.file = k;
  }
  if (a.file.empty())    { a.error  = "no input file specified"; return a; }
  if (a.formatStr.empty()) a.formatStr = pathExtname(a.file);
  try { if (!a.samplesStr.empty()) a.samples = stoi(a.samplesStr); }
  catch (...) { a.error = "\'"+a.samplesStr+"\' samples is not an integer"; return a; }
  a.format    = parseFileFormat(a.formatStr);
  a.transform = parseGraphTransform(a.transformStr);
  if (a.format   ==F::UNKNOWN) { a.error = "\'"+a.formatStr   +"\' format is not recognized";    return a; }
  if (a.transform==T::UNKNOWN) { a.error = "\'"+a.transformStr+"\' transform is not recognized"; return a; }
  if (a.samples<=0)            { a.error = "\'"+a.samplesStr+  "\' samples must be positive";    return a; }
  return a;
}




// HELP
// ----

const char* helpMessage() {
  return "For usage details, please try the following URL:\n"
  "https://github.com/puzzlef/graph-properties";
}
