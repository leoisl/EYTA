//
// Created by Leandro Ishi Soares de Lima on 17/04/17.
//

#ifndef KSGATB_STRATEGY_H
#define KSGATB_STRATEGY_H

#include "EnhanceTranscriptome.h"
#include "UnitigLinkingGraph.h"

class UnitigLinkingGraph;

namespace EnhanceTranscriptome {
    //Represents several strategies to find alternative paths. For now, we only have one...
    class Strategy {
    public:
        //return: a string with a message to the user summarizing how was the computation
        virtual void processTranscript(const MappingInfo &mappingInfo, FilteredGraph &graph,
                                       graph_t &unfilteredGraph,
                                       SubgraphVertexFilter *subgraphVertexFilter, int maxLengthAS, int maxLengthIntron,
                                       int splicingComplexityAS, int splicingComplexityIR, int k, BubbleOutputter &bubbleOutputter, const UnitigLinkingGraph &ulg,
                                       int lowerBoundTargetThreshold,
                                        //int upperBoundTargetThresholdAS, is the maximum value by default, which is (int)(mappingInfo.getMappingFW().size())
                                       int upperBoundTargetThresholdIntron) {
          throw runtime_error("Not implemented.");
        };

        //filters the graph so that only the vertices that can reach the LR nodes and that are reachable by the LR within the threshold are kept
        //this is just a graph trimming
        void getSubgraphOfInterest(FilteredGraph &graph, const MappingInfo &mappingInfo, const Vertex &artificial, int threshold,
                                   SubgraphVertexFilter *subgraphVertexFilter,
                                   int EdgeInfo::*targetWeight, int EdgeInfo::*sourceWeight);

        Strategy(int numberOfNodes, const UnitigLinkingGraph &ulg) : distances(numberOfNodes, 0), explored(numberOfNodes, false), ulg(ulg){}

        static boost::shared_ptr<Strategy> getStrategy(int strategyType, int nbOfNodes, const UnitigLinkingGraph &ulg, const map<pair<int, char>, int> &labelToIndex);

    protected:
        vector<int> distances;
        vector<bool> explored;
        const UnitigLinkingGraph &ulg;

        //all stuff concerning branches were commented out for the moment
        /*

        //get the number of branching nodes if a new vertex is added to this path
        template <class GraphType>
        int getNbBranchingNodesIfAddAVertex(const Path<GraphType> &path, const GraphType& graph) {
            if (path.nodes.size() <= 1)
                return 0;
            return path.branchingNodes + isBranching(path.nodes.back(), graph);
        }

        //Recursive function that finds all simple paths from a source s such that the path respects the constraints
        //The paths must begin in the LR, and if there is any other LR node in it, it is the last one
        //Paths are always in the forward direction, i.e. reverse graph paths are reversed
        //Do not call this function, call the getOut/InPaths()
        template <class GraphType>
        void DFSEnum(const Vertex &s, const set<Vertex> &nodesInLR, const GraphType &graph, int maxLengthOfAnAlternativeTranscript, int maxBranches,
                     list< Path<GraphType> > &paths, Path<GraphType> &currentPath, int iteration, bool isReversed, bool addOnlyPathsFinishingInTheLR) {
            //add s to the path
            explored[s]=true; //s is in the current Path

            if (iteration>0 && nodesInLR.find(s)!=nodesInLR.end()) { //if it is NOT the first iteration, and the current node is in the LR, then I found an alternative path
                //add the path
                paths.push_back(currentPath);
                if (isReversed) //reverse the path so that it represents a forward path, and not a reverse path
                    reverse(paths.back().nodes.begin(), paths.back().nodes.end());
            }

            if (iteration==0 || nodesInLR.find(s)==nodesInLR.end()) { //if iteration == 0, always explore
                //if s is not in LR, explore also
                int nbOfNeighboursExplored=0;
                for (auto outEdgesS = boost::out_edges(s, graph); outEdgesS.first != outEdgesS.second; ++outEdgesS.first) { //iterates over all out-neighbours v of s
                    auto e = *(outEdgesS.first);
                    auto v = boost::target(e, graph);

                    if (iteration == 0 && nodesInLR.find(v) != nodesInLR.end()) //if it is the first iteration, and the neighbour is in LR, skip it!
                        continue;

                    //check if v should be explored
                    int cost = (isReversed==false ? graph[e].sourceWeight : graph[e].targetWeight); //yeah, it is indeed sourceWeight - the last node of the path does not count
                    int newNbBranchingNodes = getNbBranchingNodesIfAddAVertex(currentPath, graph);
                    int newCost = (currentPath.nodes.size() > 1 ? currentPath.distance + cost : 0);
                    if (!explored[v] && //if v is not already in the path
                        newNbBranchingNodes <= maxBranches && newCost <= maxLengthOfAnAlternativeTranscript) {
                        //yes, v should be explored
                        nbOfNeighboursExplored++;

                        //Configure the currentPath accordingly
                        int oldDistance = currentPath.distance;
                        currentPath.distance = newCost;
                        int oldBranchingNodes = currentPath.branchingNodes;
                        currentPath.branchingNodes = newNbBranchingNodes;
                        currentPath.nodes.push_back(v);

                        //call DFS
                        DFSEnum(v, nodesInLR, graph, maxLengthOfAnAlternativeTranscript, maxBranches, paths, currentPath,
                                iteration + 1, isReversed, addOnlyPathsFinishingInTheLR);

                        //DesConfigure the currentPath accordingly
                        currentPath.distance = oldDistance;
                        currentPath.branchingNodes = oldBranchingNodes;
                        currentPath.nodes.pop_back();
                    }
                }
                //add this path to allPaths
                //TODO: change this to do not add the path only if ALL neighbours were explored?
                if (addOnlyPathsFinishingInTheLR==false && currentPath.nodes.size()>=2 && nbOfNeighboursExplored==0) { //if this is a subpath of another path, does not add it
                    paths.push_back(currentPath);
                    if (isReversed) //reverse the path so that it represents a forward path, and not a reverse path
                        reverse(paths.back().nodes.begin(), paths.back().nodes.end());
                }
            }

            //already explored all paths with this prefix
            explored[s]=false; //s is not in the currentPath
        }


        //add to paths all the outpaths from source in the graph respecting the constraints  maxLengthOfAnAlternativeTranscript and maxBranches
        template <class GraphType>
        void getPaths(const Vertex &source, const set<Vertex> &nodesInLR, const GraphType &graph, int maxLengthOfAnAlternativeTranscript, int maxBranches,
                      list< Path<GraphType> > &paths, bool isReversed, bool addOnlyPathsFinishingInTheLR) {
            //initialize the global variables
            fill(explored.begin(), explored.end(), false);
            Path<GraphType> currentPath(&graph);

            //find all paths
            currentPath.nodes.push_back(source);
            DFSEnum(source, nodesInLR, graph, maxLengthOfAnAlternativeTranscript, maxBranches, paths, currentPath, 0, isReversed, addOnlyPathsFinishingInTheLR);
        }

        template <class GraphType>
        void getOutPaths(const Vertex &source, const set<Vertex> &nodesInLR, GraphType &graph,
                         int maxLengthOfAnAlternativeTranscript, int maxBranches, list< Path<GraphType> > &paths, bool addOnlyPathsFinishingInTheLR) {
            getPaths(source, nodesInLR, graph, maxLengthOfAnAlternativeTranscript, maxBranches, paths, false, addOnlyPathsFinishingInTheLR);
        }
        template <class GraphType>
        void getInPaths(const Vertex &source, const set<Vertex> &nodesInLR, GraphType &graph,
                        int maxLengthOfAnAlternativeTranscript, int maxBranches, list< Path<GraphType> > &paths, bool addOnlyPathsFinishingInTheLR) {
            //get the reversed graph
            boost::reverse_graph <graph_t> reversedGraph = boost::make_reverse_graph(graph);
            getPaths(source, nodesInLR, reversedGraph, maxLengthOfAnAlternativeTranscript, maxBranches, paths, true, addOnlyPathsFinishingInTheLR);
        }
         */

        //check if both paths are vertex-disjoint
        template <class IteratorType>
        bool vertexDisjoint (const IteratorType &p1Begin, const IteratorType &p1End, const IteratorType &p2Begin, const IteratorType &p2End) {
            for (IteratorType p1It = p1Begin; p1It != p1End; ++p1It) {
                if (find(p2Begin, p2End, *p1It) != p2End)
                    return false;
            }
            return true;
        }
    };
}

#endif //KSGATB_STRATEGY_H
