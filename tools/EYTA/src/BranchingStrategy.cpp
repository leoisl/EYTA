//
// Created by Leandro Ishi Soares de Lima on 17/04/17.
//

#include "BranchingStrategy.h"
#include "Utils.h"

//commented out
//namespace EnhanceTranscriptome{
//    //TODO: maybe we should go back to the naive algorithm - the graph seems very simple here, so it is probably faster if you use a simple algorithm!!!
//    void BranchingStrategy::findAlternativePath(FilteredGraph &graph, graph_t &unfilteredGraph,
//                             SubgraphVertexFilter *subgraphVertexFilter,
//                             const Vertex &v, vector<Vertex> &alternativePath, int iteration,
//                             const set<Vertex> &targetNodes, int maxLengthPath, int branches,
//                             const MappingInfo &mappingInfo, int startIndex, bool *firstBubble,
//                             int &nbOfAlternativePaths, int maxLengthOfAnAlternativeTranscript,
//                             int maxBranches, BubbleOutputter &bubbleOutputter) {
//      /*
//       * TODO: REMOVED ALL OF THIS
//      #ifdef EYTA_DEBUG
//      //cout << "[findAlternativePath] on node " << graph[v].id << "_" << graph[v].strand << endl;
//      #endif
//      //base case of our recursive function:
//      //did we reached an alternative path to a targetNode?
//      if (iteration > 1 && targetNodes.find(v) != targetNodes.end()) {
//        //yes!
//        alternativePath.push_back(v);
//
//        bubbleOutputter.outputBubble(graph, alternativePath, mappingInfo, startIndex, firstBubble,
//                     maxLengthOfAnAlternativeTranscript - maxLengthPath);
//        nbOfAlternativePaths++;
//
//        alternativePath.pop_back();
//        return;
//      }
//
//      //add this node to the alternative path
//      alternativePath.push_back(v); //add this node to the alternative path
//      int maxDistance = maxLengthPath - graph[v].weight; //takes the weight of v out of the threshold
//
//      //check which neighbour can be added to the alternative path
//      auto adjIterator = adjacent_vertices(v, graph);
//
//      for_each(adjIterator.first, adjIterator.second, [&](const Vertex &neighbour) {
//          if (iteration == 1 && targetNodes.find(neighbour) != targetNodes.end()) {
//            //if we are in the first iteration and our neighbour is a target, ignore it
//            return;
//          }
//
//          #ifdef EYTA_DEBUG
//          //cout << "Neighbour found: " << graph[neighbour].id << "_" << graph[neighbour].strand << endl;
//          #endif
//
//          //check if the neighbour satisfies the branching condition
//          auto neighbourEdgePair = edge(v, neighbour, graph);
//          int neighbourBranches;
//          if (neighbourEdgePair.second) {
//            auto neighbourEdge = neighbourEdgePair.first;
//            neighbourBranches = branches + graph[neighbourEdge].targetIsBranching;
//            if (neighbourBranches > maxBranches) //way too many branches
//              return; //this is a return because we are in a for_each()
//          } else {
//            cout << "FATAL ERROR ON auto neighbourEdgePair = edge(v, neighbour, graph);" << endl;
//            exit(1);
//          }
//
//
//
//
//          //check if the neighbour satisfies the condition of max length of the path
//          //remove v from the graph
//          if (targetNodes.find(v) == targetNodes.end()) {
//            //will NOT enter here if v is source and in the first iteration!
//            //source node is allowed to be repeated
//            subgraphVertexFilter->remove(v);
//          }
//
//          fill(distances.begin(), distances.end(), 0);
//          int distanceToClosestLRNode = -1;
//          try {
//            //invoke dijkstra
//            dijkstra_shortest_paths(graph, neighbour, weight_map(get(&EdgeInfo::targetWeight, graph)).
//                distance_map(make_iterator_property_map(distances.begin(),
//                                                        boost::get(boost::vertex_index, graph))).
//                visitor(StopWhenReachingLRDijkstraVisitor(targetNodes, distances, maxDistance)));
//            throw LongReadWasNotReached(); //dijkstra completed, and the target was not reached - throw
//          } catch (LongReadWasReached &e) {
//            distanceToClosestLRNode = e.distance;
//          } catch (LongReadWasNotReached &e) { }
//
//          //put v back in the graph
//          if (targetNodes.find(v) == targetNodes.end()) {
//            //add v back to the graph
//            subgraphVertexFilter->add(v);
//          }
//
//          if (distanceToClosestLRNode == -1) {
//            //no path to LR nodes with the maxDistance provided - go to next neighbour
//            return; //this is a return because we are in a for_each()
//          }
//
//          //here, there is a path to LR nodes with the maxDistance provided
//          //call findAlternativePath recursively to list all paths in this subtree
//          findAlternativePath(graph, unfilteredGraph, subgraphVertexFilter, neighbour, alternativePath,
//                              iteration + 1, targetNodes, maxDistance,
//                              neighbourBranches, mappingInfo, startIndex, firstBubble, nbOfAlternativePaths,
//                              maxLengthOfAnAlternativeTranscript, maxBranches, bubbleOutputter);
//      });
//
//      //all alternative paths from this subtree were listed
//      alternativePath.pop_back();
//       */
//    }
//
//
//    void BranchingStrategy::processCluster(const vector<MappingInfo> &mappingInfo, FilteredGraph &graph, graph_t &unfilteredGraph,
//                                      SubgraphVertexFilter *subgraphVertexFilter, int maxLengthOfAnAlternativeTranscript,
//                                      int maxBranches, BubbleOutputter &bubbleOutputter, int cluster, int readId) {
//      throw runtime_error("Reimplement");
////      //TODO: take care of unmapped regions (-1 x -1)
////
////      //TODO: take care of streches!!!!
////      //TODO: streches will probably describe small variations (if any), since they are consecutive LR kmers mapping to the same SR unitig
////      //TODO: unless you have a region that is described ONLY by LR , then the stretch could be like:
////      //TODO: U1 U1 U1 U1 U1 X X X X X X X X X X X X X U1 U1 U1
////      //TODO:                <Long read unique region>
////      //TODO: we need to take care of this also!
////      //TODO: for the moment, just compact same unitig stretches, keeping the first and the last pos of each strech
////      //TODO: you use the last pos of a strech to get out of the unitig, and the first to come back when getting the constitutive kmers of a bubble
////
////      //#ifdef EYTA_DEBUG
//////          cout << "********************************************" << endl;
//////          cout << getTime() << " - Processing Mapping: ";
//////          for_each(mappingInfo.mapping.begin(), mappingInfo.mapping.end(), [&](const Vertex &v) {
//////              if (v==boost::graph_traits<graph_t>::null_vertex())
//////                cout << "X_X ";
//////              else
//////                cout << graph[v].id << "_" << graph[v].strand << " ";
//////          });
//////          cout << endl;
////      //#endif
////
////      cout << "********************************************" << endl;
////      cout << getTime() << " - Processing read: " << mappingInfo.readFileIndex << " " << mappingInfo.readIndex <<
////      endl;
////      cout << "********************************************" << endl;
////
////      //add the nodes of the LR to a set
////      set<Vertex> nodesInLR(mappingInfo.mapping.begin(), mappingInfo.mapping.end());
////      nodesInLR.erase(boost::graph_traits<graph_t>::null_vertex());
////
////      //will get the subgraph associated with this LR - all nodes that can reach and are reachable by the LR within distance D
////      //define the subgraph - the filters define the subgraphs
////      subgraphVertexFilter->reset(true); //put everyone in the graph
////
////      //connect the artificial node to all LR nodes
////      Vertex artificial = vertex(labelToIndex.at(make_pair(ARTIFICIAL_ID, 'F')), unfilteredGraph);
////      for (const auto &nodeInLR : nodesInLR) {
////        pair<Edge, bool> return_from_add_edge = add_edge(artificial, nodeInLR, unfilteredGraph);
////        if (return_from_add_edge.second) {
////          graph[return_from_add_edge.first].targetWeight = 0;
////          graph[return_from_add_edge.first].sourceWeight = 0;
////          graph[return_from_add_edge.first].sourceIsBranching = 0;
////          graph[return_from_add_edge.first].targetIsBranching = 0;
////          graph[return_from_add_edge.first].id = ARTIFICIAL_ID;
////        } else {
////          cerr << "Failed to add artificial edge!" << endl;
////          exit(1);
////        }
////        return_from_add_edge = add_edge(nodeInLR, artificial, unfilteredGraph);
////        if (return_from_add_edge.second) {
////          graph[return_from_add_edge.first].targetWeight = 0;
////          graph[return_from_add_edge.first].sourceWeight = 0;
////          graph[return_from_add_edge.first].sourceIsBranching = 0;
////          graph[return_from_add_edge.first].targetIsBranching = 0;
////          graph[return_from_add_edge.first].id = ARTIFICIAL_ID;
////        } else {
////          cerr << "Failed to add artificial edge!" << endl;
////          exit(1);
////        }
////      }
////
////      //computes the subgraph of interest
////      getSubgraphOfInterest(graph, artificial, maxBranches, subgraphVertexFilter,
////                            &EdgeInfo::targetIsBranching, &EdgeInfo::sourceIsBranching);
////      getSubgraphOfInterest(graph, artificial, maxLengthOfAnAlternativeTranscript, subgraphVertexFilter,
////                            &EdgeInfo::targetWeight, &EdgeInfo::sourceWeight);
////
////
////      //EYTA_DEBUG
////      //writeGraphToFile(graph, "graph_after_applying_constraints");
////      //
////
////      //remove the connection of the artificial to the LR nodes
////      clear_vertex(artificial, unfilteredGraph);
////
////
////      //set some vars to find the alternative paths
////      bool firstBubble[] = {true, true, true, true};
////      int i = 0;
////      int nbOfAlternativePaths = 0;
////
////
////      //find the alternative paths
////      //all arcs arriving at a LR node should have weight 0 - they do not count when building SR sequences
////      map<Edge, int> weightBackup; //to restore the previous weight
////      map<Edge, int> branchingBackup; //to restore the previous weight
////      for (const auto &nodeInLR : nodesInLR) {
////        auto inEdges = in_edges(nodeInLR, unfilteredGraph);
////        for_each(inEdges.first, inEdges.second, [&](const Edge &e) {
////            weightBackup[e] = unfilteredGraph[e].targetWeight;
////            unfilteredGraph[e].targetWeight = 0;
////            branchingBackup[e] = unfilteredGraph[e].targetIsBranching;
////            unfilteredGraph[e].targetIsBranching = 0;
////        });
////      }
////
////      for_each(mappingInfo.mapping.begin(), mappingInfo.mapping.end(), [&](const Vertex &node) {
////          BOOST_SCOPE_EXIT(&i) {
////            i++;
////          }
////          BOOST_SCOPE_EXIT_END //does not matter what happens, will increase i
////
////          if (node == boost::graph_traits<graph_t>::null_vertex()) {
////            //let's skip this node since it is invalid (unmapped)
////            return;
////          }
////
////          //TODO: what if the strech is like this:
////          //U1 U1 U1 U1 U1 X X X X X X X X X X X X X U1 U1 U1
////          //Should I go up to the last U1?
////          if (i < mappingInfo.mapping.size() - 1 && node == mappingInfo.mapping[i + 1]) {
////            //let's skip this node because it is in the middle of a strech (we just process the last node of a stretch)
////            return;
////          }
////
////          //cout << getTime() << " - Processing node " << graph[node].id << "_" << graph[node].strand << endl;
////
////          //get all the target nodes
////          set<Vertex> targetNodes(mappingInfo.mapping.begin() + i, mappingInfo.mapping.end());
////          targetNodes.erase(boost::graph_traits<graph_t>::null_vertex());
////
////          //calls the recursive function which will find the alternative paths from source to targets
////          vector<Vertex> alternativePath;
////          findAlternativePath(graph, unfilteredGraph, subgraphVertexFilter, node, alternativePath, 1,
////                              targetNodes, maxLengthOfAnAlternativeTranscript, 0,
////                              mappingInfo, i, firstBubble, nbOfAlternativePaths, maxLengthOfAnAlternativeTranscript,
////                              maxBranches, bubbleOutputter);
////      });
////
////      //restore the previous weight
////      for (const auto &nodeInLR : nodesInLR) {
////        auto inEdges = in_edges(nodeInLR, unfilteredGraph);
////        for_each(inEdges.first, inEdges.second, [&](const Edge &e) {
////            unfilteredGraph[e].targetWeight = weightBackup[e];
////            unfilteredGraph[e].targetIsBranching = branchingBackup[e];
////        });
////      }
////
////      //cout << "Found " << nbOfAlternativePaths << " alternative paths." << endl;
////
////      #ifdef EYTA_DEBUG
////      cout << "********************************************" << endl;
////      #endif
//    }
//}