//
// Created by Leandro Ishi Soares de Lima on 07/12/17.
//

#include "UnitigLinkingGraph.h"


void UnitigLinkingGraph::addEdgeCore(const UnitigLinkingGraphEdge& edge, const string &seq) {
  //if this edge already exists, we just update seqCount
  if (edges.count(edge)>0) {
    edges.find(edge)->seqCount[seq]++;
    return;
  }

  //here, we have a new edge!
  //add the edge
  edges.insert(edge);
  edges.find(edge)->seqCount[seq]++;
}

//probably more efficient in a parallel way than calling addEdge() several times
void UnitigLinkingGraph::addEdges(const vector< tuple<UnitigInfo, UnitigInfo, string> > &allULGEdgesToAddNeedsProcessing,
                                  const SolidKmerOracle &solidKmerOracle) {
  //create allULGEdgesToAdd
  vector< tuple<UnitigLinkingGraphEdge, string> > allULGEdgesToAdd;
  allULGEdgesToAdd.reserve(allULGEdgesToAddNeedsProcessing.size());
  for (const auto &triple : allULGEdgesToAddNeedsProcessing) {
    auto fromId = get<0>(triple).unitigId;
    UnitigLinkingGraphNode from(fromId, get<0>(triple).strand, solidKmerOracle.getUnitigSequence(fromId).size());
    auto toId = get<1>(triple).unitigId;
    UnitigLinkingGraphNode to(toId, get<1>(triple).strand, solidKmerOracle.getUnitigSequence(toId).size());
    UnitigLinkingGraphEdge edge(from, to);
    string seq = get<2>(triple);
    allULGEdgesToAdd.push_back(make_tuple(edge, seq));
  }


  //add all edges synchronously
  synchroAddEdge->lock();
    for (const auto &pair : allULGEdgesToAdd)
      addEdgeCore(get<0>(pair), get<1>(pair));
  synchroAddEdge->unlock();
}


//add all the hints in allHintsToAdd
void UnitigLinkingGraph::addAllHints(const vector< pair<UnitigInfo, UnitigInfo> > &allHintsToAddNeedsProcessing,
                                     const SolidKmerOracle &solidKmerOracle) {
  //create allHintsToAdd
  vector< pair<UnitigLinkingGraphNode, UnitigLinkingGraphNode> > allHintsToAdd;
  allHintsToAdd.reserve(allHintsToAddNeedsProcessing.size());
  for (const auto &pair : allHintsToAddNeedsProcessing) {
    auto fromId = pair.first.unitigId;
    UnitigLinkingGraphNode from(fromId, pair.first.strand, solidKmerOracle.getUnitigSequence(fromId).size());
    auto toId = pair.second.unitigId;
    UnitigLinkingGraphNode to(toId, pair.second.strand, solidKmerOracle.getUnitigSequence(toId).size());
    allHintsToAdd.push_back(make_pair(from, to));
  }

  synchroAddAllHints->lock();
  for (const auto &hintPair : allHintsToAdd)
    hints[hintPair.first].insert(hintPair.second);
  synchroAddAllHints->unlock();
}