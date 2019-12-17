#include "HintsDrivenEnumerationStrategy.h"
#include "Utils.h"
#include <boost/graph/breadth_first_search.hpp>

namespace EnhanceTranscriptome{
    class TooManyAlternativePaths{};


    //DEBUG - print the graph for debug
    //DEBUG - print the graph for debug
    class bfs_get_nodes_and_edges_visitor:public boost::default_bfs_visitor {
    public:
        bfs_get_nodes_and_edges_visitor(set<Vertex>& vertices):vertices(vertices){ }
        void discover_vertex(Vertex u, const FilteredGraph & g) {
          vertices.insert(u);
        }
        set<Vertex>& vertices; //vertices found so far
    };

    void printGraphForDebug(const string &prefix, FilteredGraph &graph, const Vertex &node) {
      set<Vertex> vertices; //vertices found so far
      boost::breadth_first_search(graph, node, boost::visitor(bfs_get_nodes_and_edges_visitor(vertices)));
      ofstream nodesFile;
      openFileForWriting(prefix+".nodes", nodesFile);
        for (const auto &v : vertices)
          nodesFile << GraphWriter::toStringNodeComplete(v, graph) << endl;
      nodesFile.close();

      ofstream edgesFile;
      openFileForWriting(prefix+".edges", edgesFile);
      for (const auto &from : vertices) {
        for (const auto &to : vertices) {
          auto possibleEdge = edge(from,to,graph);
          if (possibleEdge.second) {
            edgesFile << GraphWriter::toString(possibleEdge.first, graph) << endl;
          }
        }
      }
      edgesFile.close();
    }
    //DEBUG - print the graph for debug
    //DEBUG - print the graph for debug


    //TODO: maybe we should go back to the naive algorithm - the graph seems very simple here, so it is probably faster if you use a simple algorithm!!!
    void HintsDrivenEnumerationStrategy::findAlternativePath(FilteredGraph &graph, graph_t &unfilteredGraph,
                                                SubgraphVertexFilter *subgraphVertexFilter,
                                                const Vertex &v, vector<Vertex> &alternativePath, int iteration,
                                                const set<Vertex> &targetNodes, int maxDistance,
                                                const MappingInfo &mappingInfo, int startIndex, bool *firstBubble,
                                                int &nbOfAlternativePaths, int maxLengthOfAnAlternativeTranscript,
                                                int splicingComplexity, int k, BubbleOutputter &bubbleOutputter, map<Vertex, int> &hintedNodes, const string &prefix) {
      //DEBUG
      if (DEBUG_SEARCH_ALG) cout << "*******************************************************************" << endl << "[findAlternativePath] on node " << graph[v].id << "_" << graph[v].strand << endl;

      //base case of our recursive function:
      //did we reached an alternative path to a targetNode?
      if (iteration > 1 && targetNodes.find(v) != targetNodes.end()) {
        //yes!
        alternativePath.push_back(v);

        cout << "Path found between " << GraphWriter::toString(alternativePath.front(), graph) << " and " << GraphWriter::toString(alternativePath.back(), graph) << endl;

        //check if the path found is really new regarding the transcriptome or not
        if (bubbleOutputter.isAlternativePathReallyNovel(graph, alternativePath, mappingInfo)) {
          bubbleOutputter.outputBubble(graph, alternativePath, mappingInfo, startIndex, firstBubble,
                                       maxLengthOfAnAlternativeTranscript - maxDistance, prefix);
          nbOfAlternativePaths++;
        }

        if (nbOfAlternativePaths==splicingComplexity) {
          cout << "Throwing TooManyAlternativePaths..." << endl;
          throw TooManyAlternativePaths();
        }

        alternativePath.pop_back();
        return;
      }

      //2/ We flag the nodes that are hinted by s
      UnitigLinkingGraphNode ulgV(graph[v].id, graph[v].strand, graph[v].seq.size());
      if (ulg.getHints().count(ulgV)) {
        for (const auto &ulgHintedNode : ulg.getHints().at(ulgV)) {
          auto hintedNode = vertex(labelToIndex.at(make_pair(ulgHintedNode.getId(), ulgHintedNode.getStrand())), unfilteredGraph);
          hintedNodes[hintedNode]++;
        }
      }


      //check which neighbour can be added to the alternative path
      //3/ We choose ALL neighbours satisyfing the following constraints (these are the neighbours to be explored in search of an alternative path):
        //3.1/ It does not exceed the maximum length
        //3.2/ It has the highest hints
      if (DEBUG_SEARCH_ALG) cout << "Getting adjacent vertices for vertex " << v << ": " << graph[v].id << "_" << graph[v].strand << endl;
      auto adjIterator = adjacent_vertices(v, graph);
      set<Vertex> nonTargetNeighboursToBeExplored; //the non-target neighbours to be explored
      set<Vertex> targetNeighboursToBeExplored; //the target neighbours to be explored
      map<Vertex, int> maxDistanceNeighbour;
      int highestHint=-1;

      for_each(adjIterator.first, adjIterator.second, [&](const Vertex &neighbour) {
          if (DEBUG_SEARCH_ALG) cout << "Exploring neighbour: " << graph[v].id << "_" << graph[v].strand << " -> " << graph[neighbour].id << "_" << graph[neighbour].strand << endl;
          if (iteration == 1 && targetNodes.find(neighbour) != targetNodes.end()) {
            //if we are in the first iteration and our neighbour is a target, ignore it
            if (DEBUG_SEARCH_ALG) cout << "Not explored: if (iteration == 1 && targetNodes.find(neighbour) != targetNodes.end())" << endl;
            return;
          }

          //otherwise, if the neighbour is a target, then we should explore it no matter what - we found an alternative path
          if (targetNodes.find(neighbour) != targetNodes.end()) {
            targetNeighboursToBeExplored.insert(neighbour);
            if (DEBUG_SEARCH_ALG) cout << "Explored: added to the targets: if (targetNodes.find(neighbour) != targetNodes.end())" << endl;
            return;
          }

          //get the corresponding edge
          auto neighbourEdge = edge(v, neighbour, graph).first;

          //check if adding this neighbour we are still good regarding the maxDistance
          maxDistanceNeighbour[neighbour] = maxDistance - graph[neighbourEdge].targetWeight;
          if (maxDistanceNeighbour[neighbour] < 0) {
            if (DEBUG_SEARCH_ALG) cout << "Not explored: distance exceeded" << endl;
            return; //no reason to continue, we are in a path longer than expected
          }


          //now we have the checks
          //check if this neighbour is red only
          //TODO: removed this - we can consider all nodes, which here are red or purple
          /*
          if (graph[neighbour].transcriptomeCoverage>0) {
            if (DEBUG_SEARCH_ALG) cout << "Not explored: neighbour is not red" << endl;
            return; //no, the node is also blue
          }
          */



          //2. check if the shortest path going through the neighbour until a target is smaller than a threshold
          //TODO: do a reverse BFS from the targets to all neighbours here instead...
          //check if the neighbour satisfies the condition of max length of the path
          //remove v from the graph
          subgraphVertexFilter->remove(v);

          fill(distances.begin(), distances.end(), 0);
          int distanceToAnyTarget = -1;
          try {
            //invoke dijkstra
            dijkstra_shortest_paths(graph, neighbour, weight_map(get(&EdgeInfo::targetWeight, graph)).
                distance_map(make_iterator_property_map(distances.begin(),
                                                        boost::get(boost::vertex_index, graph))).
                visitor(StopWhenReachingLRDijkstraVisitor(targetNodes, distances, maxDistanceNeighbour[neighbour])));
            throw LongReadWasNotReached(); //dijkstra completed, and the target was not reached - throw
          } catch (LongReadWasReached &e) {
            distanceToAnyTarget = e.distance;
          } catch (LongReadWasNotReached &e) { }

          //put v back in the graph
          subgraphVertexFilter->add(v);


          if (distanceToAnyTarget == -1) {
            //no path to target nodes with the maxDistanceNeighbour[neighbour] provided - go to next neighbour
            if (DEBUG_SEARCH_ALG) cout << "Not explored: cant reach target" << endl;
            return; //this is a return because we are in a for_each()
          }
          //here, there is a path to the target nodes with the maxDistanceNeighbour[neighbour] provided






          //3. now, we check if this node has the highest hint
          if (hintedNodes[neighbour]>highestHint) {
            if (DEBUG_SEARCH_ALG) cout << "Explored: highest hint" << endl;
            //if the neighbour has hint larger than the highest hint, then he is the one leading the assembly now
            highestHint=hintedNodes[neighbour];
            nonTargetNeighboursToBeExplored.clear();
            nonTargetNeighboursToBeExplored.insert(neighbour);
          }else if (hintedNodes[neighbour]==highestHint) {
            if (DEBUG_SEARCH_ALG) cout << "Explored: equal hint" << endl;
            //neighbour has the same hint as the highest one
            nonTargetNeighboursToBeExplored.insert(neighbour);
          }
          else {
            //neighbour has smaller hint, so we don't care about him
            if (DEBUG_SEARCH_ALG) cout << "Not explored: hint too small" << endl;
          }
      });
      if (DEBUG_SEARCH_ALG) {
        cout << "Just after the for..." << endl;
        cout << "Highest hint = " << highestHint << endl;
        cout << "nonTargetNeighboursToBeExplored = ";
        for (const auto &node : nonTargetNeighboursToBeExplored)
          cout << graph[node].id << "_" << graph[node].strand << " ";
        cout << endl;
        cout << "targetNeighboursToBeExplored = ";
        for (const auto &node : targetNeighboursToBeExplored)
          cout << graph[node].id << "_" << graph[node].strand << " ";
        cout << endl;
      }


       //3.3/ If the highest hint happens to be 1, this means that the ancestral nodes are not contributing to the hints at all;
         //-This means that we can not connect what comes before and what comes after using a read - if we are traversing a repeat, this means that it is a repeat longer than the read size - then we are screwed;
         //-This is the hardest case in assembly - where the short reads do not give us any information at all on how to traverse a region;
           //If the region is simple, we are fine
           //If it is complex (repeat), then we are screwend
         //-What we do in this case is that we get the splicingComplexity longest neighbours - hoping that long unitigs generally represent good assemblies;
       if (highestHint==1) {
         if (DEBUG_SEARCH_ALG) cout << "Highest hint is 1, limiting the size of nonTargetNeighboursToBeExplored" << endl;

         //sort nonTargetNeighboursToBeExplored by graph[neighbour].weight in a decreasing way
         vector<Vertex> nonTargetNeighboursToBeExploredVector(nonTargetNeighboursToBeExplored.begin(), nonTargetNeighboursToBeExplored.end());
         sort(nonTargetNeighboursToBeExploredVector.begin(), nonTargetNeighboursToBeExploredVector.end(), [&](const Vertex &v1, const Vertex &v2) -> bool {
             return graph[v1].weight > graph[v2].weight;
         });

         //erase all nodes after the splicingComplexity-th
         if (nonTargetNeighboursToBeExploredVector.size()>splicingComplexity)
           nonTargetNeighboursToBeExploredVector.erase(nonTargetNeighboursToBeExploredVector.begin()+splicingComplexity, nonTargetNeighboursToBeExploredVector.end());

         //adjust nonTargetNeighboursToBeExplored
         nonTargetNeighboursToBeExplored.clear();
         nonTargetNeighboursToBeExplored.insert(nonTargetNeighboursToBeExploredVector.begin(), nonTargetNeighboursToBeExploredVector.end());
       }

      if (DEBUG_SEARCH_ALG) {
        cout << "Just after if (highestHint==1) {..." << endl;
        cout << "nonTargetNeighboursToBeExplored = ";
        for (const auto &node : nonTargetNeighboursToBeExplored)
          cout << graph[node].id << "_" << graph[node].strand << " ";
        cout << endl;
        cout << "targetNeighboursToBeExplored = ";
        for (const auto &node : targetNeighboursToBeExplored)
          cout << graph[node].id << "_" << graph[node].strand << " ";
        cout << endl;
      }


      //before calling the recursive function on the neighbour, we
      //add v to the alternative path
      alternativePath.push_back(v);
      //remove v from the graph
      subgraphVertexFilter->remove(v);



      //4/ We add the neighbour to the assembly
      //we first explore the targets, since they will compose an alternative path directly
      //call findAlternativePath recursively on all targets to list the alternative paths found
      for (const auto &neighbour : targetNeighboursToBeExplored) {
        findAlternativePath(graph, unfilteredGraph, subgraphVertexFilter, neighbour, alternativePath,
                            iteration+1, targetNodes, maxDistanceNeighbour[neighbour],
                            mappingInfo, startIndex, firstBubble, nbOfAlternativePaths,
                            maxLengthOfAnAlternativeTranscript, splicingComplexity, k, bubbleOutputter, hintedNodes, prefix);
      }

      //we now explore the non-target neighbours
      for (const auto &neighbour : nonTargetNeighboursToBeExplored) {
        findAlternativePath(graph, unfilteredGraph, subgraphVertexFilter, neighbour, alternativePath,
                            iteration+1, targetNodes, maxDistanceNeighbour[neighbour],
                            mappingInfo, startIndex, firstBubble, nbOfAlternativePaths,
                            maxLengthOfAnAlternativeTranscript, splicingComplexity, k, bubbleOutputter, hintedNodes, prefix);
      }


      //all alternative paths from this subtree were listed
      alternativePath.pop_back();
      subgraphVertexFilter->add(v);


      //unflag the nodes that are hinted by v - in order to not interfere with the other assembly
      if (ulg.getHints().count(ulgV)) {
        for (const auto &ulgHintedNode : ulg.getHints().at(ulgV)) {
          auto hintedNode = vertex(labelToIndex.at(make_pair(ulgHintedNode.getId(), ulgHintedNode.getStrand())), unfilteredGraph);
          hintedNodes[hintedNode]--;
        }
      }

      if (DEBUG_SEARCH_ALG) cout << "*******************************************************************" << endl;
    }


    void HintsDrivenEnumerationStrategy::findEvents (const MappingInfo &mappingInfo, FilteredGraph &graph,
                                                     graph_t &unfilteredGraph,
                                                     SubgraphVertexFilter *subgraphVertexFilter, int maxLengthOfAnAlternativeTranscript,
                                                     int splicingComplexity, int k, BubbleOutputter &bubbleOutputter, const UnitigLinkingGraph &ulg, const Vertex &node,
                                                     int lowerBound, int upperBound, int& nbOfAlternativePaths, int startIndex, bool *firstBubble, const string &prefix
                                                     #if DEBUG_SEARCH_ALG == 1
                                                     , ofstream& sourcesFile
                                                     #endif
    ) {
      //check if the lower bound is fine
      if (lowerBound>=mappingInfo.getMappingFW().size()) {
        #if DEBUG_SEARCH_ALG == 1
        sourcesFile << "SKIPPED: lower bound out of bounds" << endl;
        #endif
        return; // no alternative path here
      }


      //the target nodes are the nodes between the lower and the upper bound "starting" in this interval, and that are trustable and red+blue
      set<Vertex> targetNodes;
      for (int i=lowerBound; i<upperBound; i++) {
        if (mappingInfo.getMappingFW()[i] != mappingInfo.getMappingFW()[i-1]) {
          //mappingInfo.getMappingFW()[i] "starts" @ position i
          auto targetNode = mappingInfo.getMappingFW()[i];

          //check if it is trustable and red+blue
          if (unfilteredGraph[targetNode].isTrustable && unfilteredGraph[targetNode].shortReadsCoverage>0)
            targetNodes.insert(targetNode);
        }
      }

      //check if there is at least one good target
      if (targetNodes.size()==0) {
        //no need to proceed
        #if DEBUG_SEARCH_ALG == 1
        sourcesFile << "SKIPPED: no viable target" << endl;
        #endif
        return;
      }


      //all set and ready, let's try to find alternative paths between a source and the targets
      cout << getTime() << " - Processing source " << GraphWriter::toString(node, unfilteredGraph) << endl;
      #if DEBUG_SEARCH_ALG == 1
      sourcesFile << "PROCESSED" << endl;
      #endif


      //all arcs arriving at the targets should have weight 0 - they do not count when building SR sequences
      //TODO: I am just using the weight now, not the branching info
      map<Edge, int> edgeWeightBackup; //to restore the previous weight
      //map<Edge, int> branchingBackup; //to restore the previous weight
      for (const auto &targetNode :targetNodes) {
        auto inEdges = in_edges(targetNode, unfilteredGraph);
        for_each(inEdges.first, inEdges.second, [&](const Edge &e) {
            edgeWeightBackup[e] = unfilteredGraph[e].targetWeight;
            unfilteredGraph[e].targetWeight = 0;
        });
      }

      //we will also remove the weight of the source and the targets:
      map<Vertex, int> nodeWeightBackup; //to restore the previous weight
      nodeWeightBackup[node]=unfilteredGraph[node].weight;
      unfilteredGraph[node].weight=0;
      for (const auto &target : targetNodes) {
        nodeWeightBackup[target]=unfilteredGraph[target].weight;
        unfilteredGraph[target].weight=0;
      }



      //IDEA 1:
      //here, we shouldn't use the ULG - we are being very specific in the sequence bases (finding stop codons) and sequence length (multiple of 3)... the ULG do not care about keeping the sequence exact
      //do a search for the paths based on the 6 ORFs - do not include paths with stop codons
      //have a maximum local splicing complexity of 5 - you should not expect more than 5 splicing events between 2 constitutive exons
      //in this step, avoid alternative paths starting with gt and finishing with ag - we do not want to find intron retention here
      //disregard the UTRs regions of the transcript - just consider the part of the transcript between the first start codon and the last stop codon
      //-i am trying to focus on the AS evetns that have a functional impact, and these only happen in the CDS
      //if it is an AS event of a transcript, it should not change the reading frame - it must be a multiple of 3!!!
      //What I means is that the frame of the right flanking exon must be kept!!


      //IDEA 2:
      //when searching for an alternative path in the DBG in a transcript T of a gene G, we should discard all nodes that come from other genes
      //we do not want to build chimeric alternative paths that goes to other genes due to repeats - we want to keep inside our gene or gene family
      //we should also discard all nodes coming from G also, except for the source and targets of the local event we are investigating
      //we do not want to find an alternative path that is described by another transcript T' of G
      //or even in another region of T -> we are local
      //to do so, we discard all blue (or red+blue) nodes from the graph, just keeping the sources and the targets
      //but then let's say that I have a repeat R that is present in many genes. All events going through R will be lost?
      //If R is divergent enough so that we have some trustable unitigs specific to one gene, then we can recover it
      //We need an additional ULG connecting all RED trustful unitigs in a read (the REDULG)
      //This is an additional data structure to the normal ULG, connecting all trustful unitigs. The idea is that we go from a BLUETU (the source) to a REDTU using the ULG, then traverse the RED region using the REDULG and then go back to a BLUETU using the ULG
      //The reasoning is that when finding alternative paths, we are interested on the RED trustful unitigs to be used as anchors in the assembly.
      //we could have some BLUE kmers due to the edges sequences between the nodes of the REDULG, but the anchors are important, and they will be both RED and trustful
      //have a maximum local splicing complexity of 5 - you should not expect more than 5 splicing events between 2 constitutive exons
      //MCF7 graph with Homo_sapiens.GRCh38.cdna.all.fa as ref. transcriptome has:
      //26 647 224 unitigs
      //14 146 581 trustful unitigs
      //12 858 042 red-only trustful unitigs
      //10 423 352 if min_size_unitigs=10
      //35 708 125 ULG edges
      //24 631 509 red-only ULG edges


      //IDEA 3:
      //Assembly towards the longest unitigs, taking into consideration the unitigs that can be reached using read and PE-read information
      //For now, let's forget about PE-read information... Let's just use the read info
      //We now have a longest unitig linking graph (LULG)
      //For each trustful unitig, we connect it only to the top n (n=10 by default) longest unitigs it can reach
      //Note that this differs from the ULG, mainly because ULG just connects consecutive unitigs, and considers all trustful unitigs
      //Here, we look at all unitigs a source unitig s can reach, and just keep the link between s and the n longest ones (keep the representative sequence also to fill the gap!!)
      //1/ Start with a trustful unitig source s (we still have the concept of a trustful unitig):
      //2/ Scan all trustful unitigs that can be reached from s (such unitigs are the ones that are in a read together with s or in its mate-pair)
      //2.1/ If a target is among them, then we are done;
      //2.1/ Otherwise, choose the largest unitig from this set;
      //3/ Continue from point 1/, but now the source is the chosen node;
      //4/ After finding the target, you will have: s -> Unitig 1 -> Unitig 2 -> ... -> t. Now you need to fill in the gaps
      //5/ Remove all inner unitigs of the path found before searching the second path



      //IDEA 3 (simplified, to test if it could work):
      //Assembly towards the longest unitigs, taking into consideration the unitigs that can be reached using read and PE-read information
      //For now, let's forget about PE-read information... Let's just use the read info
      //We now have the FULL unitig linking graph (FULG)
      //For each trustful unitig, we connect it only to ALL unitigs it can reach
      //Note that this differs from the ULG, mainly because ULG just connects consecutive unitigs
      //Here, we look at all unitigs a source unitig s can reach using the read information (keep the representative sequence also to fill the gap!!)
      //1/ Start with a trustful unitig source s:
      //2/ Scan all trustful unitigs that can be reached from s (such unitigs are the ones that are in a read together with s or in its mate-pair)
      //2.1/ If a target is among them, then we are done;
      //2.1/ Otherwise, choose the largest unitig red-only from this set;
      //3/ Continue from point 1/, but now the source is the chosen node;
      //4/ After finding the target, you will have: s -> Unitig 1 -> Unitig 2 -> ... -> t. Now you use the arc labels to fill the gap
      //5/ Remove all inner unitigs of the path found before searching the second path
      //6/ You should search for a maximum number m of maximum paths, which is the maximum local splicing complexity

      //IDEA 4 (based on 3)
      //TODO: I SHOULD NOT CARE ABOUT THE PART OF THE TRANSCRIPTS IN THE UTRs - they do not change the CDS -> no functional AS
      //TODO: The flanking exons must be necessarily RED+BLUE, since I am building the ULG by the mapping of the short reads to the DBG...

      //Assembly towards the longest unitigs, taking into consideration the unitigs that can be reached using read and PE-read information
      //The longest unitigs do not represent repeats, and if we have reads solving repeats, we will probably have a link between the long unitigs flanking the repeat
      //For now, let's forget about PE-read information... Let's just use the read info
      //Unitig linking graph (ULG) is now like this:
      //For each trustful unitig, we connect it to the n (10 by default) largest trustful unitigs it can reach using the read information
      //We keep the representative sequence also to fill the gap!!
      //A trustful unitig is one having blue+red color with size >= 36, or if it is red only, it needs size >= 61
      //PS: blue-only nodes are currently disregarded since the flanking exons must be expressed (red) and known (blue) to find the AS event.
      //the idea is that our source and target anchors are going to be blue+red nodes, so we have to keep them
      //Plus, since blue nodes come from the transcriptome, we trust them more, thus the we are looser in its size
      //In the MCF7 dataset, we have ??? red+blue nodes
      //Red only nodes are not trustable - they come from the short reads, not from the transcriptome
      //and they are way harder to assemble due to repeats, so since we are doing an assembly towards the longest contigs, there is no bad on keeping only the 10 largest unitigs - we are going to disregard the smallest unitigs anyway
      //In the MCF7 dataset, we have 23.8M red-only nodes, and this is the distribution of their size:
      //@see https://www.dropbox.com/sh/n4t9f8nj6nr5f1a/AAAHsTapmIARHR8ibYXsjHWra?dl=0
      //after the smallest unitig sizes, we have a peak @ 61 and a peak @100 (read size)
      //that is why I am using right now 61

      //Alg:
      //1/ Start with a blue+red trustful unitig source s (the left flanking exon)
      //2/ Scan all trustful unitigs that can be reached from s (such unitigs are the ones that are in a read together with s or in its mate-pair)
      //2.1/ If a target is among them, then we found one path;
      //2.2/ Otherwise, choose the largest unitig from this set;
      //2.3/ Continue from point 1/, but now the source is the chosen node;
      //3/ After finding the target, you will have: s -> Unitig 1 -> Unitig 2 -> ... -> t. Now you need to fill in the gaps

      //4/ List all paths like this... We probably won't have many
      //4.0/ How to list all paths like this?? I will list only one for now...
      //4.1/ Use splicing complexity??

      //IDEA 5 (based on 4)
      //Each trustful unitig needs connections to Blue+Red unitigs and connections to Red-only unitigs
      //Like this, we can go from the flanking exon 1 to the Red-only path and then back to the flanking exon 2








      //IDEA 6
      //ULG is now edges + hints
      //Trustful unitigs connect ONLY to the consecutive trustful unitigs in a short read
      //If in a short read we have TU1xxxxxTU2xxxxxTU3, we will have the following 2 edges: TU1 -> TU2 e TU2 -> TU3
      //A trustful unitig u has hints: the set of unitigs S such that for each v \in S there is a short read connecting u -> v
      //One idea was to give weight to hints: if there are 50 reads connecting u -> v, then this hint should have weight 50
      //The problem is that when you enter a repeat node, this node is contained in way more reads than the nodes of the transcript you are building
      //So, it might give you a wrong hint: his hints are going to be preferred because they have very high weight, and the previous hints, given by unitigs you were assembling before,
      //will be disconsidered, but you should consider the hints old nodes more than the hint of a node that just arrived!!
      //That is why hints do not have weight

      //Algorithm for assembling (for now, just where the transcript is the exclusion isoform and we try to find an inclusion isoform):
      //1/ Start with a blue+red trustful unitig source s (the left flanking exon)
      //We assume this left flanking exon is not contained in other genes - if it is, we might (or not) assemble chimeric stuff
      //2/ We flag the nodes that are hinted by s
      //3/ We choose ALL neighbours satisyfing the following constraints (these are the neighbours to be explored in search of an alternative path):
      //3.1/ It is RED;
      //3.2/ Has the highest hint (not sure if this is good...);
      //3.3/ If the highest hint happens to be 1, this means that the ancestral nodes are not contributing to the hints at all;
      //-This means that we can not connect what comes before and what comes after using a read - if we are traversing a repeat, then we are screwed;
      //-This is the hardest case in assembly - where the short reads do not give us any information at all on how to traverse a region;
      //-What we do in this case is that we get the splicingComplexity longest neighbours - hoping that long unitigs generally represent good assemblies;
      //3.4/ TODO: add the length contraint somewhere here?
      //4/ We add the neighbour to the assembly
      //4.1/ Let's say v=neighbour
      //5/ We flag the nodes that are hinted by v
      //6/ We now choose the next node to enter in the assembly (one of v's neighbour)
      //6.1/ If a target is in v's neighbours, we put it in the path and output the path found
      //6.2/ We loop to step 3 to find all red neighbours to continue enumerating the paths

      //TODO:
      // 2. add paired-end hints to the ULG
      // 3. find other types of events (not only exclusion x inclusion isoform)
      // I am not sure if it is nice to remove all BLUE nodes...






      //calls the recursive function which will find the alternative paths from source to targets
      vector<Vertex> alternativePath;
      int nbOfAlternativePathsLocal=0;
      map<Vertex, int> hintedNodes;
      try {
        findAlternativePath(graph, unfilteredGraph, subgraphVertexFilter, node, alternativePath, 1,
                            targetNodes, maxLengthOfAnAlternativeTranscript,
                            mappingInfo, startIndex, firstBubble, nbOfAlternativePathsLocal,
                            maxLengthOfAnAlternativeTranscript,
                            splicingComplexity, k, bubbleOutputter, hintedNodes, prefix);
      }catch (const TooManyAlternativePaths& e){}
      nbOfAlternativePaths+=nbOfAlternativePathsLocal;

      //restore the previous weight
      for (const auto &targetNode : targetNodes) {
        auto inEdges = in_edges(targetNode, unfilteredGraph);
        for_each(inEdges.first, inEdges.second, [&](const Edge &e) {
            unfilteredGraph[e].targetWeight = edgeWeightBackup[e];
        });
      }

      //restore the weight of the source and the targets:
      unfilteredGraph[node].weight=nodeWeightBackup[node];
      for (const auto &target : targetNodes)
        unfilteredGraph[target].weight=nodeWeightBackup[target];
}


void HintsDrivenEnumerationStrategy::processTranscript(const MappingInfo &mappingInfo, FilteredGraph &graph,
                                             graph_t &unfilteredGraph,
                                             SubgraphVertexFilter *subgraphVertexFilter, int maxLengthAS, int maxLengthIntron,
                                             int splicingComplexityAS, int splicingComplexityIR, int k, BubbleOutputter &bubbleOutputter, const UnitigLinkingGraph &ulg,
                                             int lowerBoundTargetThreshold,
                                             //int upperBoundTargetThresholdAS, is the maximum value by default, which is (int)(mappingInfo.getMappingFW().size())
                                             int upperBoundTargetThresholdIntron) {

      //#ifdef EYTA_DEBUG
      //          cout << "********************************************" << endl;
      //          cout << getTime() << " - Processing Mapping: ";
      //          for_each(mappingInfo.mapping.begin(), mappingInfo.mapping.end(), [&](const Vertex &v) {
      //              if (v==boost::graph_traits<graph_t>::null_vertex())
      //                cout << "X_X ";
      //              else
      //                cout << graph[v].id << "_" << graph[v].strand << " ";
      //          });
      //          cout << endl;
      //#endif

      cout << "********************************************" << endl;
      cout << getTime() << " - Processing transcript: " << mappingInfo.getReadIndex() << endl;
      cout << "********************************************" << endl;
      //TODO : THIS METHOD MUST BE SYNCED
      bubbleOutputter.prepareOutput(mappingInfo);

      /*
       * TODO: skip this?
      //will get the subgraph associated with this LR - all nodes that can reach and are reachable by the LR within distance D
      //define the subgraph - the filters define the subgraphs
      subgraphVertexFilter->reset(true); //put everyone in the graph

      //connect the artificial node to all LR nodes
      Vertex artificial = vertex(labelToIndex.at(make_pair(ARTIFICIAL_ID, 'F')), unfilteredGraph);
      for (const auto &nodeInLR : mappingInfo.getUniqueNodes()) {
        pair<Edge, bool> return_from_add_edge = add_edge(artificial, nodeInLR, unfilteredGraph);
        if (return_from_add_edge.second) {
          graph[return_from_add_edge.first].targetWeight = 0;
          graph[return_from_add_edge.first].sourceWeight = 0;
          graph[return_from_add_edge.first].edgeWeight = 0;
          graph[return_from_add_edge.first].sourceIsBranching = 0;
          graph[return_from_add_edge.first].targetIsBranching = 0;
          graph[return_from_add_edge.first].label = "";
          graph[return_from_add_edge.first].id = ARTIFICIAL_ID;
        } else {
          cout << "Failed to add artificial edge!" << endl;
          exit(1);
        }
        return_from_add_edge = add_edge(nodeInLR, artificial, unfilteredGraph);
        if (return_from_add_edge.second) {
          graph[return_from_add_edge.first].targetWeight = 0;
          graph[return_from_add_edge.first].sourceWeight = 0;
          graph[return_from_add_edge.first].edgeWeight = 0;
          graph[return_from_add_edge.first].sourceIsBranching = 0;
          graph[return_from_add_edge.first].targetIsBranching = 0;
          graph[return_from_add_edge.first].label = "";
          graph[return_from_add_edge.first].id = ARTIFICIAL_ID;
        } else {
          cout << "Failed to add artificial edge!" << endl;
          exit(1);
        }
      }

      //computes the subgraph of interest - based on the branches
      //getSubgraphOfInterest(graph, mappingInfo, artificial, maxBranches, subgraphVertexFilter,
      //                      &EdgeInfo::targetIsBranching, &EdgeInfo::sourceIsBranching);

      //computes the subgraph of interest - based on the weights
      //TODO: I should fix the weights to include the edge weight
      //getSubgraphOfInterest(graph, mappingInfo, artificial, maxLengthOfAnAlternativeTranscript, subgraphVertexFilter,
      //                      &EdgeInfo::targetWeight, &EdgeInfo::sourceWeight);

      //EYTA_DEBUG
      //writeGraphToFile(graph, "graph_after_applying_constraints");
      //

      //remove the connection of the artificial to the LR nodes
      clear_vertex(artificial, unfilteredGraph);

             */


      //set some vars to find the alternative paths
      bool firstBubble[] = {true, true, true, true};
      int i = 0;
      int nbOfAlternativePaths = 0;


      //find the alternative paths
      //TODO: also work on getMappingRC()?


      //DEBUG
      //write all sources in fasta format to file
      #if DEBUG_SEARCH_ALG == 1
      ofstream sourcesFile;
      openFileForWriting("read_"+to_string(mappingInfo.getReadIndex())+"_sources.fa", sourcesFile);
      #endif

      for_each(mappingInfo.getMappingFW().begin(), mappingInfo.getMappingFW().end(), [&](const Vertex &node) {
          BOOST_SCOPE_EXIT(&i) {
            i++;
          }
          BOOST_SCOPE_EXIT_END //does not matter what happens, will increase i

          //do some initial checks to verify if we should proceed

          //just process the last node of a stretch
          if (i < mappingInfo.getMappingFW().size() - 1 && node == mappingInfo.getMappingFW()[i + 1]) {
            //let's skip this node because it is in the middle of a strech (we just process the last node of a stretch)
            return;
          }

          //DEBUG
          //print the graph for debug so we can see it!
          //printGraphForDebug("read_"+to_string(mappingInfo.getReadFileIndex())+"_"+to_string(mappingInfo.getReadIndex())+"_graph_source_"+GraphWriter::toString(node, graph), graph, node);
          #if DEBUG_SEARCH_ALG == 1
          sourcesFile << ">" << GraphWriter::toString(node, unfilteredGraph) << endl << graph[node].seq << endl;
          #endif

          //check if the start is trustable
          if (unfilteredGraph[node].isTrustable==false) {
            //lets skip it
            #if DEBUG_SEARCH_ALG == 1
            sourcesFile << "SKIPPED: not trustable" << endl;
            #endif
            return;
          }

          //check if the node is also expressed - if it is red
          //the flanking exons should be red and blue
          if (unfilteredGraph[node].shortReadsCoverage==0) {
            //not red - lets skip it
            #if DEBUG_SEARCH_ALG == 1
            sourcesFile << "SKIPPED: not red" << endl;
            #endif
            return;
          }


          //now it seems we should proceed - at least the constraints for the starting node is fine
          //check the other constraints and find AS events
          findEvents(mappingInfo, graph, unfilteredGraph, subgraphVertexFilter, maxLengthAS, splicingComplexityAS, k, bubbleOutputter, ulg, node,
                     i+lowerBoundTargetThreshold, (int)(mappingInfo.getMappingFW().size()), nbOfAlternativePaths, i, firstBubble, "Novel_AS"
          #if DEBUG_SEARCH_ALG == 1
          , sourcesFile
          #endif
          );

          //check the other constraints and find intron events
          //basically, for introns, we try to find only ONE path between the flanking exons (way too much repeats... it is hard to find events here...)
            //however, it is not a random path, it is the path that makes more sense possible (regarding the hints / unitig lengths)
          //the length is a lot larger, but the targets are a lot smaller
          findEvents(mappingInfo, graph, unfilteredGraph, subgraphVertexFilter, maxLengthIntron, splicingComplexityIR, k, bubbleOutputter, ulg, node,
                     i+lowerBoundTargetThreshold, min((int)(mappingInfo.getMappingFW().size()), i+upperBoundTargetThresholdIntron), nbOfAlternativePaths, i, firstBubble, "Novel_IR"
          #if DEBUG_SEARCH_ALG == 1
          , sourcesFile
          #endif
          );

     });
      #if DEBUG_SEARCH_ALG == 1
      sourcesFile.close();
      #endif
      cout << "Found " << nbOfAlternativePaths << " alternative paths." << endl;

      //cleanup everything related to the mapping of this transcript
      //TODO : THIS METHOD MUST BE SYNCED
      bubbleOutputter.finalizeOutput(mappingInfo);

      #ifdef EYTA_DEBUG
      cout << "********************************************" << endl;
      #endif
    }
}