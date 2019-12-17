//
// Created by Leandro Ishi Soares de Lima on 17/04/17.
//

#include "BranchingPlusBubbleGeneratorStrategy.h"
#include "Utils.h"
#include "EnhanceTranscriptomeDefs.h"
#include <set>

using namespace std;

//commented out
//namespace EnhanceTranscriptome {
//    void BranchingPlusBubbleGeneratorStrategy::processCluster(const vector<MappingInfo> &mappingInfo, FilteredGraph &doNotUseThisGraph, graph_t &graph,
//                   SubgraphVertexFilter *subgraphVertexFilter, int maxLengthOfAnAlternativeTranscript,
//                   int maxBranches, BubbleOutputter &bubbleOutputter, int cluster, int readId) {
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
////      cout << "********************************************" << endl;
////      cout << getTime() << " - Processing read: " << mappingInfo.readFileIndex << " " << mappingInfo.readIndex <<
////      endl;
////      cout << "********************************************" << endl;
////
////      //add the nodes of the LR to a set
////      set<Vertex> nodesInLR(mappingInfo.mapping.begin(), mappingInfo.mapping.end());
////      nodesInLR.erase(boost::graph_traits<graph_t>::null_vertex());
////
////      //now, we get all the out and in-paths from and to a node of the LR that respects the constraints
////      map<Vertex,list<Path> > outPaths;
////      map<Vertex,list<Path> > inPaths;
////      for (const auto &nodeInLR : nodesInLR) {
////        getOutPaths(nodeInLR, nodesInLR, graph, maxLengthOfAnAlternativeTranscript, maxBranches, outPaths[nodeInLR], false);
////        getInPaths(nodeInLR, nodesInLR, graph, maxLengthOfAnAlternativeTranscript, maxBranches, inPaths[nodeInLR], false);
////      }
////
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////      //EYTA_DEBUG: PRINT THE PATHS THAT WERE FOUND
////      /*
////      cout << "-------------------------------------------------------" << endl;
////      cout << "EYTA_DEBUG: PRINT THE PATHS THAT WERE FOUND" << endl;
////      cout << "outPaths:" << endl;
////      for (const auto &nodeInLR : nodesInLR) {
////        cout << "-------------------------------------------------------" << endl;
////        cout << "Source node = " << GraphWriter::toString(nodeInLR, graph) << endl;
////        for (const auto &path : outPaths[nodeInLR])
////          cout << GraphWriter::toString(path, graph) << endl;
////        cout << "-------------------------------------------------------" << endl;
////      }
////      cout << "inPaths:" << endl;
////      for (const auto &nodeInLR : nodesInLR) {
////        cout << "-------------------------------------------------------" << endl;
////        cout << "Source node = " << GraphWriter::toString(nodeInLR, graph) << endl;
////        for (const auto &path : inPaths[nodeInLR])
////          cout << GraphWriter::toString(path, graph) << endl;
////        cout << "-------------------------------------------------------" << endl;
////      }
////      cout << "-------------------------------------------------------" << endl;
////       */
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////      //////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////      //set some vars to find the alternative paths
////      bool firstBubble[] = {true, true, true, true};
////      int i = 0;
////      int nbOfAlternativeTrivialPaths = 0;
////      int nbOfAlternativeRepeatContainingPaths = 0;
////      set< vector<Vertex> > allOutputPaths; //will store all the output paths so that we do not output repeated paths
////
////
////      //now, list all paths such that they start and end in the LR - these are called the trivial bubbles
////      //i.e. the ones not going through a complex region (they were found by the algorithm with less than maxBranches branches and respect the length constraints)
////      //as we output these trivial bubbles, we remove them also because they will not be needed
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
////          //the core of the listing bubbles is here
////          {
////            list<Path> &outPathsOfThisNode = outPaths[node];
////            for (auto outPathsOfThisNodeIt = outPathsOfThisNode.begin(); outPathsOfThisNodeIt!=outPathsOfThisNode.end();) {
////              if (nodesInLR.find((*outPathsOfThisNodeIt).nodes.front())!=nodesInLR.end() && nodesInLR.find((*outPathsOfThisNodeIt).nodes.back())!=nodesInLR.end()) {
////                //this is a trivial bubble - output and remove
////                //output
////                if (bubbleOutputter.outputBubble(graph, (*outPathsOfThisNodeIt).nodes, mappingInfo, i, firstBubble,
////                                                  (*outPathsOfThisNodeIt).distance, &allOutputPaths))
////                  nbOfAlternativeTrivialPaths++;
////
////                //remove
////                outPathsOfThisNodeIt=outPathsOfThisNode.erase(outPathsOfThisNodeIt);
////              }else {
////                //increment by hand
////                ++outPathsOfThisNodeIt;
////              }
////            }
////          }
////          //do the same with the in-paths - just for the sake of removing them
////          {
////            list<Path> &inPathsOfThisNode = inPaths[node];
////            for (auto inPathsOfThisNodeIt = inPathsOfThisNode.begin(); inPathsOfThisNodeIt!=inPathsOfThisNode.end();) {
////              if (nodesInLR.find((*inPathsOfThisNodeIt).nodes.front())!=nodesInLR.end() && nodesInLR.find((*inPathsOfThisNodeIt).nodes.back())!=nodesInLR.end()) {
////                //this is a trivial bubble - remove
////                inPathsOfThisNodeIt=inPathsOfThisNode.erase(inPathsOfThisNodeIt);
////              }else {
////                //increment by hand
////                ++inPathsOfThisNodeIt;
////              }
////            }
////          }
////      });
////
////      //now let's try to list the trivial paths formed by a concatenation of an outpath and an inpath
////      //this only happens when the last node of the outpath is exactly the same as the first one of the inpath, and the rest of the nodes are vertex-disjoint
////      i=0;
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
////          //for each pair <node, targetNode>, list the complicated paths
////          //TODO: THIS IS WAY, WAY TOO SLOW... WHAT TO DO?
////          for (const auto &targetNode : targetNodes) {
////            //for each outPath of node and each inPath of the target node
////            for (const auto &outPath : outPaths[node]) {
////              for (const auto &inPath : inPaths[targetNode]) {
////                //check if the paths are vertex disjoint, apart from their extremal nodes
////                if (vertexDisjoint(outPath.nodes.begin(), outPath.nodes.end() - 1, inPath.nodes.begin() + 1,
////                                   inPath.nodes.end()) == false)
////                  //both paths are not vertex-disjoint, they can not represet an alternative path
////                  continue;
////
////                //check if the last node of the outpath is equal to the first node of the inpath
////                if (outPath.nodes.back() != inPath.nodes.front())
////                  //not equal, I should not list them here
////                  continue;
////
////                //output
////                vector<Vertex> bothParts(outPath.nodes);
////                for_each(inPath.nodes.begin() + 1, inPath.nodes.end(), [&](const Vertex &node) {
////                    bothParts.push_back(node);
////                });
////                if (bubbleOutputter.outputBubble(graph, bothParts, mappingInfo, i, firstBubble,
////                                                                      outPath.distance+inPath.distance-graph[outPath.nodes.back()].weight,
////                                                                      &allOutputPaths))
////                  nbOfAlternativeTrivialPaths++;
////              }
////            }
////          }
////      });
////
////      //now let's try to list the complicated bubbles, that might contain repeats
////      i=0;
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
////          //for each pair <node, targetNode>, list the complicated paths
////          for (const auto& targetNode : targetNodes) {
////            //for each outPath of node and each inPath of the target node
////            for(const auto& outPath : outPaths[node]) {
////              for(const auto& inPath : inPaths[targetNode]) {
////                //check if the paths are vertex disjoint, apart from their extremal nodes
////                if (vertexDisjoint(outPath.nodes.begin(), outPath.nodes.end()-1, inPath.nodes.begin()+1, inPath.nodes.end())==false)
////                  //both paths are not vertex-disjoint, they can not represet an alternative path
////                  continue;
////
////                //check if the last node of the outpath is different of the first node of the inpath
////                if (outPath.nodes.back()==inPath.nodes.front())
////                  //they are equal, I should not list them here
////                  continue;
////
////
////                //create the graph without the nodes from the two paths
////                //except the last node of the outPath and the first node of the inPath
////                auto subgraphEdgeFilter = new SubgraphEdgeFilterRemoveSet(&graph);
////                auto subgraphVertexFilter = new SubgraphVertexFilterRemoveSet(&graph, subgraphEdgeFilter);
////                FilteredGraphRemoveSet graphWithoutTheTwoPaths(graph, SubgraphEdgeFilterForFG<SubgraphEdgeFilterRemoveSet>(subgraphEdgeFilter),
////                                            SubgraphVertexFilterForFG<SubgraphVertexFilterRemoveSet>(subgraphVertexFilter));
////                for (const auto &node : outPath.nodes)
////                  subgraphVertexFilter->remove(node);
////                for (const auto &node : inPath.nodes)
////                  subgraphVertexFilter->remove(node);
////                subgraphVertexFilter->add(outPath.nodes.back());
////                subgraphVertexFilter->add(inPath.nodes.front());
////
////
////                //in the given graph, finds the shortest path from the outPath to the inPath
////                //TODO: is the shortest path the best option here?
////                fill(distances.begin(), distances.end(), 0); //TODO: change this for a map?
////                int distanceToInPath = -1;
////                set<Vertex> targetNodes;
////                targetNodes.insert(inPath.nodes.front());
////                try {
////                  //invoke dijkstra
////                  dijkstra_shortest_paths(graphWithoutTheTwoPaths, outPath.nodes.back(), weight_map(get(&EdgeInfo::targetWeight, graphWithoutTheTwoPaths)).
////                      distance_map(make_iterator_property_map(distances.begin(),
////                                                              boost::get(boost::vertex_index, graphWithoutTheTwoPaths))).
////                      //TODO: is the size maxLengthOfAnAlternativeTranscript fine here?
////                      visitor(StopWhenReachingLRDijkstraVisitor(targetNodes, distances, maxLengthOfAnAlternativeTranscript)));
////                  throw LongReadWasNotReached(); //dijkstra completed, and the target was not reached - throw
////                } catch (LongReadWasReached &e) {
////                  distanceToInPath = e.distance;
////                } catch (LongReadWasNotReached &e) { }
////
////                if (distanceToInPath == -1) {
////                  //no path from the outPath to the inPath was found
////                  continue;
////                }
////
////
////                //here we have a path that we have to output, which is composed of:
////                //OutPath + (N * distanceToInPath) + InPath
////                if (bubbleOutputter.outputRepeatContainingBubble(graph, outPath.nodes, distanceToInPath, inPath.nodes,
////                                                                                  mappingInfo, i, firstBubble,
////                                                                                      outPath.distance+distanceToInPath+inPath.distance,
////                                                                                      &allOutputPaths))
////                  nbOfAlternativeRepeatContainingPaths++;
////              }
////            }
////          }
////      });
////
////      //cout << "Found " << nbOfAlternativeTrivialPaths << " trivial alternative paths." << endl;
////      //cout << "Found " << nbOfAlternativeRepeatContainingPaths << " repeat containing alternative paths." << endl;
////
////      #ifdef EYTA_DEBUG
////      cout << "********************************************" << endl;
////      #endif
//    }
//}