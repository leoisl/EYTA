//
// Created by Leandro Ishi Soares de Lima on 07/12/17.
//

#ifndef EYTA_HintsDrivenEnumerationStrategy_H
#define EYTA_HintsDrivenEnumerationStrategy_H

#include "Strategy.h"
#include "UnitigLinkingGraph.h"

#define DEBUG_SEARCH_ALG 0

namespace EnhanceTranscriptome {
    class HintsDrivenEnumerationStrategy : public Strategy {
    private:
        const map<pair<int, char>, int> &labelToIndex;

        void findEvents (const MappingInfo &mappingInfo, FilteredGraph &graph,
                         graph_t &unfilteredGraph,
                         SubgraphVertexFilter *subgraphVertexFilter, int maxLengthOfAnAlternativeTranscript,
                         int splicingComplexity, int k, BubbleOutputter &bubbleOutputter, const UnitigLinkingGraph &ulg, const Vertex &node,
                         int lowerBound, int upperBound, int& nbOfAlternativePaths, int startIndex, bool *firstBubble, const string &prefix
        #if DEBUG_SEARCH_ALG == 1
            , ofstream& sourcesFile
        #endif
        );

        //TODO: maybe we should go back to the naive algorithm - the graph seems very simple here, so it is probably faster if you use a simple algorithm!!!
        void findAlternativePath(FilteredGraph &graph, graph_t &unfilteredGraph,
                                                                 SubgraphVertexFilter *subgraphVertexFilter,
                                                                 const Vertex &v, vector<Vertex> &alternativePath, int iteration,
                                                                 const set<Vertex> &targetNodes, int maxDistance,
                                                                 const MappingInfo &mappingInfo, int startIndex, bool *firstBubble,
                                                                 int &nbOfAlternativePaths, int maxLengthOfAnAlternativeTranscript,
                                                                 int splicingComplexity, int k, BubbleOutputter &bubbleOutputter, map<Vertex, int> &hintedNodes, const string &prefix);
    public:
        virtual void processTranscript(const MappingInfo &mappingInfo, FilteredGraph &graph,
                                       graph_t &unfilteredGraph,
                                       SubgraphVertexFilter *subgraphVertexFilter, int maxLengthAS, int maxLengthIntron,
                                       int splicingComplexityAS, int splicingComplexityIR, int k, BubbleOutputter &bubbleOutputter, const UnitigLinkingGraph &ulg,
                                       int lowerBoundTargetThreshold,
                                        //int upperBoundTargetThresholdAS, is the maximum value by default, which is (int)(mappingInfo.getMappingFW().size())
                                       int upperBoundTargetThresholdIntron);

        HintsDrivenEnumerationStrategy(int numberOfNodes, const UnitigLinkingGraph &ulg, const map<pair<int, char>, int> &labelToIndex) : Strategy(numberOfNodes, ulg),
                                                                                              labelToIndex(labelToIndex) { }
    };

}


#endif //EYTA_HintsDrivenEnumerationStrategy_H
