#pragma once
#include <string>
#include <string_view>
#include <unordered_map>
#include <iostream>
#include <stdexcept>
#include "inc/transpose.hxx"

using std::string;
using std::string_view;
using std::to_string;
using std::unordered_map;
using std::random_device;
using std::runtime_error;
using std::mt19937_64;



#pragma region OPTIONS
/**
 * @brief
 * The dynamic graphs are generated in batch updates of specified size B,
 * either in absolute terms, or in terms of the number of edges of the
 * original graphs (like 0.001|E|). The generator should accept as parameter
 * the mix of edge deletions and insertions, mix of vertex deletions and
 * insertions (and at what rate the vertices increase/scale), whether
 * duplicate edge deletions/insertions appear in the batch update, the
 * nature of the batch update (uniformly random, preferential attachment,
 * planted partition model, ..., match the nature of the graph), and the
 * set of constraints to be satisfied (limit min/max degree, maintain degree
 * distribution, maintain community structures, k-core, diameter, ...).
 * Further, in addition to generating batch updates from the base static
 * graphs (single-batch update), it should be able to generate multi-batch
 * updates, where each batch update comes one after the other (in a sequence)
 * to form a sequence of graphs that represent the dynamic graph.
 */


/** Command-line options for the program. */
struct Options {
  unordered_map<string, string> params;
  vector<string> transforms;
};

/**
 * Check if a string is a boolean.
 * @param x string to check
 * @returns true if the string is a boolean, false otherwise
 */
inline bool isBool(string_view x) {
  return x=="0" || x=="1" || x=="true" || x=="false";
}


/**
 * Check if a string is an integer.
 * @param x string to check
 * @returns true if the string is an integer, false otherwise
 */
inline bool isInt(string_view x) {
  char *end;
  strtol(x.data(), &end, 10);
  return end == x.data() + x.size();
}


/**
 * Check if a string is a double.
 * @param x string to check
 * @returns true if the string is a double, false otherwise
 */
inline bool isDouble(string_view x) {
  char *end;
  strtod(x.data(), &end);
  return end == x.data() + x.size();
}


/**
 * Extract the extension from a path.
 * @param path path to extract extension from
 * @returns extension of the path
 */
inline string_view pathExtname(string_view path) {
  auto   idx = path.rfind('.');
  return idx==string::npos? "" : path.substr(idx);
}


/**
 * Read the command line options, and fill in default values.
 * @param argc number of arguments
 * @param argv array of arguments
 * @returns the options
 */
inline Options readOptions(int argc, char **argv) {
  Options o;
  for (int i=1; i<argc; ++i) {
    string k = argv[i];
    if (k=="--help") o.params["help"] = "1";
    else if (k=="--input-graph")     o.params["input-graph"]     = argv[++i];
    else if (k=="--input-format")    o.params["input-format"]    = argv[++i];
    else if (k=="--input-transform"){ while (i+1<argc && argv[i+1][0]!='-') o.transforms.push_back(argv[++i]);}
    else if (k=="--output-dir")      o.params["output-dir"]    = argv[++i];
    else if (k=="--output-prefix")   o.params["output-prefix"] = argv[++i];
    else if (k=="--output-format")   o.params["output-format"] = argv[++i];
    else if (k=="--batch-size")       o.params["batch-size"]       = argv[++i];
    else if (k=="--batch-size-ratio") o.params["batch-size-ratio"] = argv[++i];
    else if (k=="--edge-insertions")  o.params["edge-insertions"]  = argv[++i];
    else if (k=="--edge-deletions")   o.params["edge-deletions"]   = argv[++i];
    else if (k=="--allow-duplicate-edges") o.params["allow-duplicate-edges"] = "1";
    else if (k=="--vertex-insertions")     o.params["vertex-insertions"]  = argv[++i];
    else if (k=="--vertex-deletions")      o.params["vertex-deletions"]   = argv[++i];
    else if (k=="--vertex-growth-rate")    o.params["vertex-growth-rate"] = argv[++i];
    else if (k=="--allow-duplicate-vertices") o.params["allow-duplicate-vertices"] = "1";
    else if (k=="--update-nature") o.params["update-nature"] = argv[++i];
    else if (k=="--probability-distribution") o.params["probability-distribution"] = argv[++i];
    else if (k=="--min-degree")    o.params["min-degree"]    = argv[++i];
    else if (k=="--max-degree")    o.params["max-degree"]    = argv[++i];
    else if (k=="--max-diameter")  o.params["max-diameter"]  = argv[++i];
    else if (k=="--preserve-degree-distribution") o.params["preserve-degree-distribution"] = "1";
    else if (k=="--preserve-communities")         o.params["preserve-communities"] = "1";
    else if (k=="--preserve-k-core")              o.params["preserve-k-core"] = argv[++i];
    else if (k=="--multi-batch") o.params["multi-batch"] = argv[++i];
    else if (k=="--seed") o.params["seed"] = argv[++i];
  }
  return o;
}
#pragma endregion




#pragma region HELP
/**
 * Maybe the user needs help.
 * @returns something helpful
 */
inline const char* helpMessage() {
  // Input formats: edgelist,matrix-market,snap-temporal
  // Input transforms: transpose,unsymmetrize,symmetrize,loop-deadends,loop-vertices,clear-weights,set-weights
  // Output formats: edgelist
  const char *message =
  "Usage: graph-generate [OPTIONS]\n"
  "\n"
  "Options:\n"
  "  --input-graph <file>           Path to the input static graph file.\n"
  "  --input-format <format>        Format of the input static graph file.\n"
  "  --input-transform <transforms> Transformations to apply to the input graph.\n"
  "  --output-dir <directory>       Directory to save the generated dynamic graphs.\n"
  "  --output-prefix <prefix>       Prefix for the generated dynamic graph files.\n"
  "  --output-format <format>       Format of the generated batch updates.\n"
  "\n"
  "Batch Size:\n"
  "  --batch-size <size>           Absolute size of each batch update.\n"
  "  --batch-size-ratio <ratio>    Size of each batch update as a fraction of the total edges (e.g., 0.001).\n"
  "\n"
  "Edge Operations:\n"
  "  --edge-insertions <percentage>   Percentage of edge insertions in each batch.\n"
  "  --edge-deletions <percentage>    Percentage of edge deletions in each batch.\n"
  "  --allow-duplicate-edges          Allow duplicate edge operations within a batch.\n"
  "\n"
  "Vertex Operations:\n"
  "  --vertex-insertions <percentage> Percentage of vertex insertions in each batch.\n"
  "  --vertex-deletions <percentage>  Percentage of vertex deletions in each batch.\n"
  "  --vertex-growth-rate <rate>      Rate at which the number of vertices increases.\n"
  "  --allow-duplicate-vertices       Allow duplicate vertex operations within a batch.\n"
  "\n"
  "Batch Update Nature:\n"
  "  --update-nature <nature>         Nature of the batch updates. Options:\n"
  "                                   uniform: Uniformly random updates.\n"
  "                                   preferential: Preferential attachment model.\n"
  "                                   planted: Planted partition model.\n"
  "                                   match: Match the nature of the original graph.\n"
  " --probability-distribution <f(x)>  Probability distribution function for the batch updates.\n"
  "\n"
  "Constraints:\n"
  "  --min-degree <degree>            Minimum degree constraint for the graph.\n"
  "  --max-degree <degree>            Maximum degree constraint for the graph.\n"
  "  --max-diameter <diameter>        Ensure the diameter of the graph does not exceed the specified value.\n"
  "  --preserve-degree-distribution   Ensure the degree distribution is maintained.\n"
  "  --preserve-communities           Preserve community structures.\n"
  "  --preserve-k-core <k>            Ensure the graph maintains a k-core structure.\n"
  "\n"
  "Multi-Batch Updates:\n"
  "  --multi-batch <num>              Number of contiguous batch updates to generate.\n"
  "\n"
  "Miscellaneous:\n"
  "  --seed <seed>                    Seed for random number generator (for reproducibility).\n"
  "  --help                           Display this help and exit.\n"
  "\n";
  return message;
}
#pragma endregion
