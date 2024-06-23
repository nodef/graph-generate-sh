#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include "inc/main.hxx"
#include "options.hxx"
#include "inc/properties.hxx"

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
#ifdef OPENMP
void handleInputFormat(const string& inputFormat, DiGraph<int, int, int>& graph, const string& inputGraph) {
  if (inputFormat == "matrix-market") {
    readMtxOmpW(graph, inputGraph.c_str());
  } else if (inputFormat == "edgelist") {
    // handle edgelist format
  } else if (inputFormat == "snap-temporal"){
    // handle snap-temporal format
  } else {
    throw runtime_error("Unknown input format: " + inputFormat);
  }
}
#else
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
#endif

/**
* @brief Handle the input transformation (transpose,unsymmetrize,symmetrize,loop-deadends,loop-vertices,clear-weights,set-weights) for the graph.
* @param inputTransform The input transformation to apply.
* @param graph The graph object to be transformed.
* @throws runtime_error if the input transformation is unknown.
*/

#ifdef OPENMP
void handleInputTransform(const string& inputTransform, DiGraph<int, int, int>& graph) {
  if (inputTransform == "");
  else if (inputTransform == "transpose") {
    graph = transposeOmp(graph);
  } else if (inputTransform == "symmetrize") {
    graph = symmetrizeOmp(graph);
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
#else
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
#endif

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
* @param probabilityDistribution The probability distribution function to use for the update.
* @param updateNature The update nature to apply.
* @param graph The graph object to be updated.
* @param rng The random number generator.
* @param batchSize The size of the batch update.
* @param edgeDeletions The fraction of edges to be deleted.
* @param edgeInsertions The fraction of edges to be inserted.
* @param weights The weights according to the given distribution must be filled in here.
* @param allowDuplicateEdges Allow duplicate edges in the batch update.
* @throws runtime_error if the update nature is unknown.
*/
void handleUpdateNature(const string& probabilityDistribution, const string& updateNature, DiGraph<int, int, int>& graph, mt19937_64& rng, size_t batchSize, double edgeDeletions, double edgeInsertions, vector<double>&weights, bool allowDuplicateEdges = true) {
  vector<tuple<int, int, int>> insertions, deletions;
  if (updateNature == "") {
    customUpdate(probabilityDistribution ,rng, graph, batchSize, edgeInsertions, edgeDeletions, insertions, deletions, allowDuplicateEdges,weights);
  } else if (updateNature == "uniform") {
    uniformUpdate(rng, graph, batchSize, edgeInsertions, edgeDeletions, insertions, deletions, allowDuplicateEdges,weights);
  } else if (updateNature == "preferential") {
    preferentialUpdate(rng, graph, batchSize, edgeInsertions, edgeDeletions, insertions, deletions, allowDuplicateEdges,weights);
  } else if (updateNature == "planted") {
    // Handle planted update 
  } else if (updateNature == "match") {
    // Handle match update 
  } else {
    throw runtime_error("Unknown update nature: " + updateNature);
  }
  applyBatchUpdateU(graph, deletions, insertions);
}

template <class G>
void writeGraphPropertiesToJSON(const G& graph, const string& filename,double divergence_list) {
  auto order = graph.order();
  auto size = graph.size();  
  auto graphDensity = ::density(graph);
  auto degrees = degreeDistribution<DiGraph<int, int, int>, int>(graph);
  double min_degree = numeric_limits<double>::max();
  double max_degree = 0.0;
  double sum_degrees = 0.0;
  for (auto degree : degrees) {
    min_degree = min(min_degree, static_cast<double>(degree));
    max_degree = max(max_degree, static_cast<double>(degree));
    sum_degrees += degree;
  }
  auto avg_degree = sum_degrees / order;
  auto scc = tarjanSCC(graph);
  auto diameter = getDiameter(graph);
  ofstream file(filename);
  file << "{" << endl;
  file << "    \"order\": " << order << "," << endl;
  file << "    \"size\": " << size << "," << endl;
  file << "    \"density\": " << fixed << setprecision(5) << graphDensity << "," << endl;
  file << "    \"degree\": {" << endl;
  file << "        \"min\": " << min_degree << "," << endl;
  file << "        \"max\": " << max_degree << "," << endl;
  file << "        \"avg\": " << fixed << setprecision(5) << avg_degree << endl;
  file << "    }," << endl;
  file << "    \"diameter\": " << diameter << "," << endl;
  if(divergence_list!=1e9)
  file << "    \"KLD\": " << divergence_list << "," << endl;
  file << "    \"degreeDistribution\": [";
  for (size_t i = 0; i < degrees.size(); ++i) {
    file << degrees[i];
    if (i != degrees.size() - 1) {
      file << ", ";
    }
  }
  file << "]," << endl;
  file << "    \"scc\": " << scc << endl;
  file << "}" << endl;
}
#pragma endregion

#pragma region MAIN HANDLER
/**
 * @brief Handles the processing of options passed to the program.
 * @param options A map containing the options and their corresponding values.
 */
void handleOptions(const Options& options) {
  auto startTime = timeNow();
  if (options.params.count("help")) {
    cout << helpMessage();
    return;
  }
  vector<string> inputTransform = options.transforms;
  string inputGraph = options.params.count("input-graph") ? options.params.at("input-graph") : "";
  string inputFormat = options.params.count("input-format") ? options.params.at("input-format") : "";
  string outputDir = options.params.count("output-dir") ? options.params.at("output-dir") : "";
  string outputPrefix = options.params.count("output-prefix") ? options.params.at("output-prefix") : "";
  string outputFormat = options.params.count("output-format") ? options.params.at("output-format") : string("edgelist");
  string propertiesFile = options.params.count("show-properties") ? options.params.at("show-properties") : "";
  int64_t batchSize = options.params.count("batch-size") ? stoll(options.params.at("batch-size")) : 0;
  double batchSizeRatio = options.params.count("batch-size-ratio") ? stod(options.params.at("batch-size-ratio")) : 0.0;
  double edgeInsertions = options.params.count("edge-insertions") ? stod(options.params.at("edge-insertions")) : 0.0;
  double edgeDeletions = options.params.count("edge-deletions") ? stod(options.params.at("edge-deletions")) : 0.0;
  bool allowDuplicateEdges = options.params.count("allow-duplicate-edges");
  double vertexInsertions = options.params.count("vertex-insertions") ? stod(options.params.at("vertex-insertions")) : 0.0;
  double vertexDeletions = options.params.count("vertex-deletions") ? stod(options.params.at("vertex-deletions")) : 0.0;
  double vertexGrowthRate = options.params.count("vertex-growth-rate") ? stod(options.params.at("vertex-growth-rate")) : 0.0;
  bool allowDuplicateVertices = options.params.count("allow-duplicate-vertices");
  string probabilityDistribution = options.params.count("probability-distribution") ? options.params.at("probability-distribution") : "";
  string updateNature = options.params.count("update-nature") ? options.params.at("update-nature") : "";
  int64_t minDegree = options.params.count("min-degree") ? stoll(options.params.at("min-degree")) : 0;
  int64_t maxDegree = options.params.count("max-degree") ? stoll(options.params.at("max-degree")) : 0;
  int64_t maxDiameter = options.params.count("max-diameter") ? stoll(options.params.at("max-diameter")) : 0;
  bool preserveDegreeDistribution = options.params.count("preserve-degree-distribution");
  bool preserveCommunities = options.params.count("preserve-communities");
  int64_t preserveKCore = options.params.count("preserve-k-core") ? stoll(options.params.at("preserve-k-core")) : 0;
  int64_t multiBatch = options.params.count("multi-batch") ? stoll(options.params.at("multi-batch")) : 1;
  random_device rd;
  int64_t seed = options.params.count("seed") ? stoll(options.params.at("seed")) : rd();
  DiGraph<int, int, int> graph;
  int counter = 0;
  ofstream outputFile;
  mt19937_64 rng(seed);
  checkInputFile(inputGraph);
  handleInputFormat(inputFormat, graph, inputGraph);
  printf("Read graph: %.3f seconds\n", duration(startTime) / 1000.0);
  if(propertiesFile != "") writeGraphPropertiesToJSON(graph, propertiesFile + outputPrefix + "_" + to_string(0),1e9);
  for(int i=0; i<inputTransform.size(); i++) {
    handleInputTransform(inputTransform[i], graph);
    printf("Perform transform %s: %.3f seconds\n", inputTransform[i].c_str(), duration(startTime) / 1000.0);
  }
  while (multiBatch--) {
    if (batchSize == 0) batchSize = graph.size() * batchSizeRatio;
    vector< double> weights;
    handleUpdateNature(probabilityDistribution, updateNature, graph, rng, batchSize, edgeDeletions, edgeInsertions,weights, allowDuplicateEdges);
    std::vector<double> normalised_weights_actual = normalize(weights);
    std::map<size_t, size_t> inDegreeDistribution;
    calculateInDegreeDistribution<DiGraph<int, int, int>, int>(graph, inDegreeDistribution);
    std::vector<double> normalised_weights_real= degreeDistributionToProbability(inDegreeDistribution);
    printf("Perform batch update %d: %.3f seconds\n", counter+1, duration(startTime) / 1000.0);
    createOutputFile(outputDir, outputPrefix, ++counter, outputFile);
    writeOutput(outputFile, graph);
    printf("Write batch update %d: %.3f seconds\n", counter, duration(startTime) / 1000.0);
    double divergence=0;
    try {
        divergence = KLDivergence(normalised_weights_real, normalised_weights_actual);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
   if(propertiesFile != "") writeGraphPropertiesToJSON(graph, propertiesFile + outputPrefix + "_" + to_string(counter),divergence);

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
