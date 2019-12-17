#ifndef EYTA_TESTS_H
#define EYTA_TESTS_H

#include <gatb/gatb_core.hpp>
#include <string>
#include "global.h"
#include "Utils.h"
#include "SolidKmerOracle.h"
#include <set>
using namespace std;

//void printAllOutgoingEdgeToFile (const string &filename);
class SolidKmerOracle;
class EYTATester{
public:
    //check is the relative error removal procedure really works
    //check if the file produced by Graph::relativeErrorRemoval() (graph.removed_edges) reflects the edges removed in the graph
    static void relativeErrorRemovalTest(Graph &graph, int nbCores, const string &prefix, double relativeCutoff);
    static void checkIfAddingBackTranscriptomicEdgeIsCorrect(const set<string> &SRGraphRemovedEdges,
                                                             const set<string> &edgesInTheTranscriptomicGraph,
                                                             const set<string> &correctlyRemovedEdges);
    static void checkRemoveEdgesFromFile(Graph &graph, const string &prefixSR, const string &prefixTr, int nbCores, const set<string> &correctlyRemovedEdges);

    //check if the solid kmer oracle is fine
    static void testSolidKmerOracle(SolidKmerOracle* solidKmerOracle);

};


#endif