//
// Created by Leandro Ishi Soares de Lima on 17/04/17.
//

#include "Strategy.h"
#include "global.h"
#include "BranchingStrategy.h"
#include "BranchingPlusBubbleGeneratorStrategy.h"
#include "BranchingNaiveStrategy.h"
#include "HintsDrivenEnumerationStrategy.h"

namespace EnhanceTranscriptome {
    boost::shared_ptr<Strategy> Strategy::getStrategy(int strategyType, int nbOfNodes, const UnitigLinkingGraph &ulg, const map<pair<int, char>, int> &labelToIndex) {
      return boost::shared_ptr<Strategy>{new HintsDrivenEnumerationStrategy(nbOfNodes, ulg, labelToIndex)};
      //TODO: commented out to compile
      /*
      switch (strategyType) {
        case BRANCHING_NAIVE_STRATEGY:
          return boost::shared_ptr<Strategy>{new BranchingNaiveStrategy(nbOfNodes)};
        case BRANCHING_STRATEGY:
          return boost::shared_ptr<Strategy>{new BranchingStrategy(nbOfNodes, labelToIndex)};
        case BRANCHING_PLUS_BG_STRATEGY:
          return boost::shared_ptr<Strategy>{new BranchingPlusBubbleGeneratorStrategy(nbOfNodes)};
        default:
          throw runtime_error("Strategy not defined");
      }
       */
    }


    //filters the graph so that only the vertices that can reach the LR nodes and that are reachable by the LR within the threshold are kept
    //this is just a graph trimming
    void Strategy::getSubgraphOfInterest(FilteredGraph &graph, const MappingInfo &mappingInfo, const Vertex &artificial, int threshold,
                                                  SubgraphVertexFilter *subgraphVertexFilter,
                                                  int EdgeInfo::*targetWeight, int EdgeInfo::*sourceWeight) {
      //invoke dijkstra
      //TODO: check this
      fill(distances.begin(), distances.end(), 0);
      set <Vertex> verticesReachableByTheSourceWithinMaxDistance;
      try {
        dijkstra_shortest_paths(graph, artificial, weight_map(get(targetWeight, graph)).
            distance_map(make_iterator_property_map(distances.begin(),
                                                    boost::get(boost::vertex_index, graph))).
            visitor(StopWhenVeryDistantFromSourceDijkstraVisitor(distances, threshold,
                                                                 verticesReachableByTheSourceWithinMaxDistance)));
      } catch (TooDistant &e) { }


      //EYTA_DEBUG
      /*
      cout << "verticesReachableByTheSourceWithinMaxDistance = {";
      for(const auto &v : verticesReachableByTheSourceWithinMaxDistance) {
        cout << graph[v].id << "_" << graph[v].strand << ", ";
      }
      cout << "}" << endl;
      */
      //EYTA_DEBUG



      //get the reversed graph
      boost::reverse_graph <FilteredGraph> reversedGraph = boost::make_reverse_graph(graph);

      //EYTA_DEBUG
      /*
      for(int i=0;i<num_vertices(graph);i++) {
        Vertex v = vertex(i, unfilteredGraph);
        cout << "Outedges of " << toString(v, unfilteredGraph) << ":" << endl;
        auto outEdges = out_edges(v, unfilteredGraph);
        for_each(outEdges.first, outEdges.second, [&](const Edge &e) {
            cout << toString(source(e, unfilteredGraph), unfilteredGraph) << " " << toString(target(e, unfilteredGraph), unfilteredGraph) << " " <<
                unfilteredGraph[e].sourceWeight << " " << unfilteredGraph[e].targetWeight <<endl;
        });
      }
      */

      //invoke the reverse dijkstra
      fill(distances.begin(), distances.end(), 0);
      set <Vertex> verticesThatCanReachTargetWithinMaxDistance;
      try {
        boost::dijkstra_shortest_paths(reversedGraph, artificial,
                                       boost::weight_map(boost::get(sourceWeight, reversedGraph)).
                                           distance_map(boost::make_iterator_property_map(distances.begin(),
                                                                                          boost::get(
                                                                                              boost::vertex_index,
                                                                                              reversedGraph))).
                                           visitor(StopWhenVeryDistantFromSourceDijkstraVisitor(distances, threshold,
                                                                                                verticesThatCanReachTargetWithinMaxDistance)));
      } catch (TooDistant &e) { }

      //EYTA_DEBUG
      /*
      cout << "verticesThatCanReachTargetWithinMaxDistance = {";
      for(const auto &v : verticesThatCanReachTargetWithinMaxDistance) {
        cout << graph[v].id << "_" << graph[v].strand << ", ";
      }
      cout << "}" << endl;
       */
      //EYTA_DEBUG

      //remove artificial from the sets - it is an artificial node
      verticesReachableByTheSourceWithinMaxDistance.erase(artificial);
      verticesThatCanReachTargetWithinMaxDistance.erase(artificial);

      //the final subgraph is an intersection of all marked vertices
      //computes the intersection
      vector <Vertex> intersection(verticesReachableByTheSourceWithinMaxDistance.size());
      vector<Vertex>::iterator itIntersection;
      itIntersection = set_intersection(verticesReachableByTheSourceWithinMaxDistance.begin(),
                                        verticesReachableByTheSourceWithinMaxDistance.end(),
                                        verticesThatCanReachTargetWithinMaxDistance.begin(),
                                        verticesThatCanReachTargetWithinMaxDistance.end(), intersection.begin());
      intersection.resize(itIntersection - intersection.begin());


      //computes the final subgraph
      subgraphVertexFilter->reset(false); //remove everyone from the graph
      //add the nodes in the intersection only
      for (const auto &v : intersection)
        subgraphVertexFilter->add(v);

      //TODO: I should probably debug here by printing the subgraph to a file and checking if it is fine

      //print some info to know how good this is
      #ifdef EYTA_DEBUG
      cout << "Subgraph for this LR:" << endl;
      cout << "Graph size: " << num_vertices(graph) << endl;
      cout << "Intersection size: " << intersection.size() << endl;

      //print the unitigs ordered by size for debugging
      /*
      vector<string> allUnitigs;
      for (const auto &v : intersection)
        allUnitigs.push_back(graph[v].name);
      sort(allUnitigs.begin(), allUnitigs.end(), [](const string &s1, const string &s2) {
          return s1.length()>s2.length();
      });
      stringstream ss;
      ss << "longest_unitigs_" << mappingInfo.getReadFileIndex() << "_" << mappingInfo.getReadIndex();
      ofstream longestUnitigsFile(ss.str());
      int i=0;
      for (const auto &s : allUnitigs) {
        longestUnitigsFile << ">longest_unitig_" << i << endl;
        longestUnitigsFile << s << endl;
        i++;
      }
      longestUnitigsFile.close();
       */

      #endif
    }
}