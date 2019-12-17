//
// Created by Leandro Ishi Soares de Lima on 24/04/17.
//

#include "BranchingNaiveStrategy.h"

//commented out
//namespace EnhanceTranscriptome {
//    void BranchingNaiveStrategy::processCluster(const vector<MappingInfo> &mappingInfoCluster, FilteredGraph &graph, graph_t &unfilteredGraph,
//                        SubgraphVertexFilter *subgraphVertexFilter, int maxLengthOfAnAlternativeTranscript,
//                        int maxBranches, BubbleOutputter &bubbleOutputter, int cluster, int readId) {
//      //TODO: take care of unmapped regions (-1 x -1)
//
//      //TODO: take care of streches!!!!
//      //TODO: streches will probably describe small variations (if any), since they are consecutive LR kmers mapping to the same SR unitig
//      //TODO: unless you have a region that is described ONLY by LR , then the stretch could be like:
//      //TODO: U1 U1 U1 U1 U1 X X X X X X X X X X X X X U1 U1 U1
//      //TODO:                <Long read unique region>
//      //TODO: we need to take care of this also!
//      //TODO: for the moment, just compact same unitig stretches, keeping the first and the last pos of each strech
//      //TODO: you use the last pos of a strech to get out of the unitig, and the first to come back when getting the constitutive kmers of a bubble
//
//      //#ifdef EYTA_DEBUG
//      //          cout << "********************************************" << endl;
//      //          cout << getTime() << " - Processing Mapping: ";
//      //          for_each(mappingInfoCluster.mapping.begin(), mappingInfoCluster.mapping.end(), [&](const Vertex &v) {
//      //              if (v==boost::graph_traits<graph_t>::null_vertex())
//      //                cout << "X_X ";
//      //              else
//      //                cout << graph[v].id << "_" << graph[v].strand << " ";
//      //          });
//      //          cout << endl;
//      //#endif
//      /*
//      cout << "********************************************" << endl;
//      cout << getTime() << " - Processing read: " << mappingInfoCluster.readFileIndex << " " << mappingInfoCluster.readIndex <<
//      endl;
//      cout << "********************************************" << endl;
//       */
//
//      //add the nodes of the LR to a set
//      //mark all LR nodes
//      set<Vertex> nodesInLRsOfTheCluster;
//      for (const auto &mappingInfo : mappingInfoCluster)
//        nodesInLRsOfTheCluster.insert(mappingInfo.getUniqueNodes().begin(), mappingInfo.getUniqueNodes().end());
//
//      //now, we get all the out and in-paths from and to a node of the LR that respects the constraints
//      //set some vars to find the alternative paths
//      bool firstBubble[] = {true, true, true, true};
//      int i = 0;
//      int nbOfAlternativePaths = 0;
//
//
//      //Find all alternative paths
//      list<Path<FilteredGraph> > allAlternativePaths;
//      //For each LR node:
//      for (const auto &node : nodesInLRsOfTheCluster) {
//        //Find alternative paths by leaving the LR node and going to SR ONLY nodes, until you go back to a LR node;
//        //get the outpaths
//        list< Path<FilteredGraph> > outPathsOfThisNode;
//        getOutPaths(node, nodesInLRsOfTheCluster, graph, maxLengthOfAnAlternativeTranscript, maxBranches,
//                    outPathsOfThisNode, true);
//
//        //add the paths found to a list
//        for (const auto &path : outPathsOfThisNode) {
//          //TODO: check this - I am now just outputting the path related to SR only nodes - I am excluding the first and last nodes which also belong to LRs!
//          //TODO: check this - I am now just outputting the path related to SR only nodes - I am excluding the first and last nodes which also belong to LRs!
//          //allAlternativePaths.push_back(path);
//          Path<FilteredGraph> pathToBeAdded(path);
//          pathToBeAdded.nodes = vector<Vertex>(path.nodes.begin() + 1, path.nodes.end() - 1);
//          allAlternativePaths.push_back(pathToBeAdded);
//        }
//      }
//
//      //output all the reads and the alternative paths
//      bubbleOutputter.outputPathFromTheCluster(graph, allAlternativePaths, mappingInfoCluster, cluster, firstBubble);
//    }
//}