//
// Created by Leandro Ishi Soares de Lima on 27/01/17.
//

#ifndef KSGATB_ENHANCETRANSCRIPTOME_H
#define KSGATB_ENHANCETRANSCRIPTOME_H

#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/scope_exit.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "BoostGraphFilters.h"
#include "EnhanceTranscriptomeDefs.h"
#include "Utils.h"
#include "BubbleOutputter.h"
#include "UnitigLinkingGraph.h"

using namespace std;

namespace EnhanceTranscriptome {
    class EnhanceTranscriptome {
    private:
        map<pair<int, char>, int> labelToIndex; //label is an entry in the .nodes file and index is the index of that node in the Boost Graph. E.g.: <14, 'R'> and 500
        map<int, pair<int, char> > indexToLabel;
        graph_t graph;
        int k;
        int maxLengthAS;
        int maxLengthIntron;
        int artificialIndex;
        int splicingComplexityAS;
        int splicingComplexityIR;
    public:
        //constructor
        EnhanceTranscriptome(const string &prefix,
                             const string &outputPrefix,
                             int k,
                             int maxLengthAS,
                             int maxLengthIntron,
                             double editDistance,
                             int splicingComplexityAS,
                             int splicingComplexityIR,
                             bool outputContext,
                             int nbCores);
    };
}

#endif //KSGATB_ENHANCETRANSCRIPTOME_H
