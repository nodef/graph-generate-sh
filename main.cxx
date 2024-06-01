#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include "inc/main.hxx"
#include "options.hxx"

using namespace std;


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
 * Write a graph in the edge list format to an output file.
 * @tparam K The vertex ID type.
 * @tparam V The vertex data type.
 * @tparam E The edge weight type.
 * @param outputFile The output file stream.
 * @param graph The directed graph to write.
 * @param weighted A flag indicating whether to print edge weights.
 */
template <class K, class V, class E>
inline void writeEdgeList(ofstream& outputFile, const DiGraph<K, V, E>& graph, bool weighted=true) {
  outputFile << graph.order() << " " << graph.size() << "\n";
  graph.forEachVertex([&](K u, V d) {
    graph.forEachEdge(u, [&](K v, E w) {
      outputFile << to_string(u) << " " << to_string(v);
      if (weighted) outputFile << " " << to_string(w);
      outputFile << "\n";
    });
  });
}

/**
* @brief Write the graph to the output file.
* @param outputFile The ofstream object for the output file.
* @param graph The graph object to be written.
*/
void writeOutput(ofstream& outputFile, const DiGraph<int, int, int>& graph) {
  writeEdgeList(outputFile, graph);
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


#pragma region MAIN
/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int main(int argc, char **argv) {
  Options o = readOptions(argc, argv);
  handleOptions(o);
  return 0;
}
#pragma endregion
