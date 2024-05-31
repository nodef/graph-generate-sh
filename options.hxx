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
typedef unordered_map<string, string> Options;


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
    if (k=="--help") o["help"] = "1";
    else if (k=="--input-graph")     o["input-graph"]     = argv[++i];
    else if (k=="--input-format")    o["input-format"]    = argv[++i];
    else if (k=="--input-transform") o["input-transform"] = argv[++i];
    else if (k=="--output-dir")      o["output-dir"]    = argv[++i];
    else if (k=="--output-prefix")   o["output-prefix"] = argv[++i];
    else if (k=="--output-format")   o["output-format"] = argv[++i];
    else if (k=="--batch-size")       o["batch-size"]       = argv[++i];
    else if (k=="--batch-size-ratio") o["batch-size-ratio"] = argv[++i];
    else if (k=="--edge-insertions")  o["edge-insertions"]  = argv[++i];
    else if (k=="--edge-deletions")   o["edge-deletions"]   = argv[++i];
    else if (k=="--allow-duplicate-edges") o["allow-duplicate-edges"] = "1";
    else if (k=="--vertex-insertions")     o["vertex-insertions"]  = argv[++i];
    else if (k=="--vertex-deletions")      o["vertex-deletions"]   = argv[++i];
    else if (k=="--vertex-growth-rate")    o["vertex-growth-rate"] = argv[++i];
    else if (k=="--allow-duplicate-vertices") o["allow-duplicate-vertices"] = "1";
    else if (k=="--update-nature") o["update-nature"] = argv[++i];
    else if (k=="--min-degree")    o["min-degree"]    = argv[++i];
    else if (k=="--max-degree")    o["max-degree"]    = argv[++i];
    else if (k=="--max-diameter")  o["max-diameter"]  = argv[++i];
    else if (k=="--preserve-degree-distribution") o["preserve-degree-distribution"] = "1";
    else if (k=="--preserve-communities")         o["preserve-communities"] = "1";
    else if (k=="--preserve-k-core")              o["preserve-k-core"] = argv[++i];
    else if (k=="--multi-batch") o["multi-batch"] = argv[++i];
    else if (k=="--seed") o["seed"] = argv[++i];
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

#pragma region HANDLE OPTIONS

#pragma region SUB HANDLERS

/**
* @brief Check if the input graph file exists and can be opened.
* @param inputGraph The path to the input graph file.
* @throws runtime_error if the input graph file cannot be found.
*/
void checkInputFile(const string& inputGraph) {
  ifstream inputFile;
  inputFile.open(inputGraph, std::ios::in); 
  if (!inputFile) {
    throw runtime_error("Input graph file not found: " + inputGraph);
  }
}

/**
* @brief Handle the input format for reading the graph.
* @param inputFormat The input format (edgelist, matrix-market, snap-temporal).
* @param graph The graph object to be populated.
* @param inputGraph The path to the input graph file.
* @throws runtime_error if the input format is unknown.
*/
void handleInputFormat(const string& inputFormat, DiGraph<int, int, int>& graph, const string& inputGraph) {
  if (inputFormat == "matrix-market") {
    readMtxW(graph, inputGraph.c_str());
  } else if (inputFormat == "edgelist") {
    // handle edgelist format
  } else if (inputFormat == "snap-temporal"){
    // handle snap-temporal format
  } else {
    throw runtime_error("Unknown input format: " + inputFormat);
  }
}

/**
* @brief Handle the input transformation (transpose,unsymmetrize,symmetrize,loop-deadends,loop-vertices,clear-weights,set-weights) for the graph.
* @param inputTransform The input transformation to apply.
* @param graph The graph object to be transformed.
* @throws runtime_error if the input transformation is unknown.
*/
void handleInputTransform(const string& inputTransform, DiGraph<int, int, int>& graph) {
  if (inputTransform == "");
  else if (inputTransform == "transpose") {
    graph = transpose(graph);
  } else if (inputTransform == "symmetrize") {
    graph = symmetrize(graph);
  } else if (inputTransform == "unsymmetrize") {
    // graph = unsymmetrize(graph);
  } else if (inputTransform == "loop-deadends") {
    // graph = loopDeadends(graph);
  } else if (inputTransform == "loop-vertices") {
    // graph = loopVertices(graph);
  } else if (inputTransform == "clear-weights") {
    // graph = clearWeights(graph);
  } else if (inputTransform == "set-weights") {
    // graph = setWeights(graph);
  } else {
    throw runtime_error("Unknown input transform: " + inputTransform);
  }
}

/**
* @brief Create an output file with a specified prefix and counter.
* @param outputDir The directory path for the output file.
* @param outputPrefix The prefix for the output file name.
* @param counter The counter value for the output file name.
* @param outputFile The ofstream object for the output file.
* @throws runtime_error if the output file cannot be created.
*/
void createOutputFile(const string& outputDir, const string& outputPrefix, int& counter, ofstream& outputFile) {
  string outputFileName = outputDir + outputPrefix + "_" + to_string(counter);
  outputFile.open(outputFileName);
  if (!outputFile) {
    throw runtime_error("Failed to create file: " + outputFileName);
  }
}

/**
* @brief Write the graph to the output file.
* @param outputFile The ofstream object for the output file.
* @param graph The graph object to be written.
*/
void writeOutput(ofstream& outputFile, const DiGraph<int, int, int>& graph) {
  write(outputFile, graph, true);
  outputFile.close();
}

/**
* @brief Handle the update nature (uniform, preferential, planted, match) for batch updates.
* @param updateNature The update nature to apply.
* @param graph The graph object to be updated.
* @param rng The random number generator.
* @param batchSize The size of the batch update.
* @param edgeDeletions The fraction of edges to be deleted.
* @param edgeInsertions The fraction of edges to be inserted.
* @param allowDuplicateEdges Allow duplicate edges in the batch update.
* @throws runtime_error if the update nature is unknown.
*/
void handleUpdateNature(const string& updateNature, DiGraph<int, int, int>& graph, mt19937_64& rng, size_t batchSize, double edgeDeletions, double edgeInsertions, bool allowDuplicateEdges = true) {
  vector<tuple<int, int, int>> insertions, deletions;
  if (updateNature == "uniform") {
    uniformUpdate(rng, graph, batchSize, edgeInsertions, edgeDeletions, insertions, deletions, allowDuplicateEdges);
  } else if (updateNature == "preferential") {
    preferentialUpdate(rng, graph, batchSize, edgeInsertions, edgeDeletions, insertions, deletions, allowDuplicateEdges);
  } else if (updateNature == "planted") {
    // Handle planted update 
  } else if (updateNature == "match") {
    // Handle match update 
  } else {
    throw runtime_error("Unknown update nature: " + updateNature);
  }
  applyBatchUpdateU(graph, deletions, insertions);
}
#pragma endregion

#pragma region MAIN HANDLER
/**
 * @brief Handles the processing of options passed to the program.
 * @param options A map containing the options and their corresponding values.
 */
void handleOptions(const Options& options) {
  if (options.count("help")) {
    cout << helpMessage();
    return;
  }
  string inputGraph = options.count("input-graph") ? options.at("input-graph") : "";
  string inputFormat = options.count("input-format") ? options.at("input-format") : "";
  string inputTransform = options.count("input-transform") ? options.at("input-transform") : "";
  string outputDir = options.count("output-dir") ? options.at("output-dir") : "";
  string outputPrefix = options.count("output-prefix") ? options.at("output-prefix") : "";
  string outputFormat = options.count("output-format") ? options.at("output-format") : string("edgelist");
  int64_t batchSize = options.count("batch-size") ? stoll(options.at("batch-size")) : 0;
  double batchSizeRatio = options.count("batch-size-ratio") ? stod(options.at("batch-size-ratio")) : 0.0;
  double edgeInsertions = options.count("edge-insertions") ? stod(options.at("edge-insertions")) : 0.0;
  double edgeDeletions = options.count("edge-deletions") ? stod(options.at("edge-deletions")) : 0.0;
  bool allowDuplicateEdges = options.count("allow-duplicate-edges");
  double vertexInsertions = options.count("vertex-insertions") ? stod(options.at("vertex-insertions")) : 0.0;
  double vertexDeletions = options.count("vertex-deletions") ? stod(options.at("vertex-deletions")) : 0.0;
  double vertexGrowthRate = options.count("vertex-growth-rate") ? stod(options.at("vertex-growth-rate")) : 0.0;
  bool allowDuplicateVertices = options.count("allow-duplicate-vertices");
  string updateNature = options.count("update-nature") ? options.at("update-nature") : string("uniform");
  int64_t minDegree = options.count("min-degree") ? stoll(options.at("min-degree")) : 0;
  int64_t maxDegree = options.count("max-degree") ? stoll(options.at("max-degree")) : 0;
  int64_t maxDiameter = options.count("max-diameter") ? stoll(options.at("max-diameter")) : 0;
  bool preserveDegreeDistribution = options.count("preserve-degree-distribution");
  bool preserveCommunities = options.count("preserve-communities");
  int64_t preserveKCore = options.count("preserve-k-core") ? stoll(options.at("preserve-k-core")) : 0;
  int64_t multiBatch = options.count("multi-batch") ? stoll(options.at("multi-batch")) : 1;
  random_device rd;
  int64_t seed = options.count("seed") ? stoll(options.at("seed")) : rd();

  DiGraph<int, int, int> graph;
  checkInputFile(inputGraph);
  handleInputFormat(inputFormat, graph, inputGraph);
  handleInputTransform(inputTransform, graph);
  int counter = 0;
  ofstream outputFile;
  createOutputFile(outputDir, outputPrefix, counter, outputFile);
  writeOutput(outputFile, graph);
  mt19937_64 rng(seed);
  while (multiBatch--) {
    if (batchSize == 0) batchSize = graph.size() * batchSizeRatio;
    handleUpdateNature(updateNature, graph, rng, batchSize, edgeDeletions, edgeInsertions, allowDuplicateEdges);
    createOutputFile(outputDir, outputPrefix, ++counter, outputFile);
    writeOutput(outputFile, graph);
  }
}
#pragma endregion
#pragma endregion