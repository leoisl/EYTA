//
// Created by Leandro Ishi Soares de Lima on 07/12/17.
//

#ifndef EYTA_UNITIGLINKINGGRAPH_H
#define EYTA_UNITIGLINKINGGRAPH_H

#include <string>
#include <sstream>
#include <algorithm>
#include <set>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <vector>
#include <tuple>
#include "SolidKmerOracle.h"
/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/


using namespace std;


//some foward declarations
class UnitigLinkingGraph;
class SolidKmerOracle;
struct UnitigInfo;

class UnitigLinkingGraphNode {
public:
    //constructor
    UnitigLinkingGraphNode(int id, char strand, int size):id(id), strand(strand), size(size){}

    //constructor for serialization only
    UnitigLinkingGraphNode(){}

    //an ULG node has an id and a strand: nodes with the same ids, but different strands are still different
    bool operator< (const UnitigLinkingGraphNode &that) const {
      if (this->id != that.id)
        return this->id < that.id;
      return this->strand < that.strand;
    }
    //an ULG node has an id and a strand: nodes with the same ids, but different strands are still different
    bool operator!= (const UnitigLinkingGraphNode &that) const {
      return (this->id!=that.id) || (this->strand!=that.strand);
    }

    //getters
    int getId() const { return id; }
    char getStrand() const { return strand; }
    int getSize() const { return size; }

    //toString
    string toString() const {
      stringstream ss;
      ss << id << " " << strand << " " << size;
      return ss.str();
    }

    //serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & id;
      ar & strand;
      ar & size;
    }

private:
    int id;
    char strand;
    int size; //unitig size
    friend class UnitigLinkingGraph;
};
BOOST_CLASS_VERSION(UnitigLinkingGraphNode, 1)



class UnitigLinkingGraph;
class UnitigLinkingGraphEdge {
public:
    UnitigLinkingGraphEdge(const UnitigLinkingGraphNode& from, const UnitigLinkingGraphNode&to) : from(from), to(to) {}

    //constructor for serialization only
    UnitigLinkingGraphEdge(){}

    //getters
    const UnitigLinkingGraphNode& getFrom() const { return from; }
    const UnitigLinkingGraphNode& getTo() const { return to; }

    //get the sequence with the highest count
    string getRepresentativeSequence() const {
      return
          max_element(seqCount.begin(), seqCount.end(), [&](const pair<string, int> &p1, const pair<string, int> &p2) {
              return p1.second < p2.second;
          })->first;
    }

    //add this seq to this edge. If the seq does not exist, create it with value = 1
    //if it already exists, just increments
    void addSeq(const string &seq) {
      seqCount[seq]++;
    }

    //comparator for edge
    bool operator< (const UnitigLinkingGraphEdge &that) const {
      if (this->from != that.from)
        return this->from < that.from;
      return this->to < that.to;
    }

    //toString
    string toString() const {
      stringstream ss;
      ss << from.toString() << " " << getRepresentativeSequence() << " " << to.toString();
      return ss.str();
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & from;
      ar & to;
      ar & seqCount;
    }

private:
    UnitigLinkingGraphNode from;
    UnitigLinkingGraphNode to;
    mutable map<string, int> seqCount;
    friend class UnitigLinkingGraph;
};
BOOST_CLASS_VERSION(UnitigLinkingGraphEdge, 1)


//Represents an abstraction of the cDBG storing as nodes the trustful unitigs
//and as edges the links between the nodes that are validated by reads
//we also have the hints given by reads/PE reads
class UnitigLinkingGraph {
public:
    UnitigLinkingGraph():synchroAddEdge(System::thread().newSynchronizer()), synchroAddAllHints(System::thread().newSynchronizer()){}

    //getters
    const set<UnitigLinkingGraphEdge>& getEdges() const {
      return edges;
    }
    const map<UnitigLinkingGraphNode, set<UnitigLinkingGraphNode> >& getHints() const {
      return hints;
    };

    //add edges in batch: probably more efficient in a parallel way than calling addEdge() several times
    void addEdges(const vector< tuple<UnitigInfo, UnitigInfo, string> > &allULGEdgesToAdd, const SolidKmerOracle &solidKmerOracle);

    //add hints in batch
    void addAllHints(const vector< pair<UnitigInfo, UnitigInfo> > &allHintsToAddNeedsProcessing,
                     const SolidKmerOracle &solidKmerOracle);

    //toString
    string toStringEdges() const {
      stringstream ss;
      for (const auto &edge : edges)
        ss << edge.toString() << endl;
      return ss.str();
    }

    //toString hints
    string toStringHints() const {
      stringstream ss;
      for (const auto &hintPair : hints) {
        ss << hintPair.first.toString() << " -> ";
        for (const auto &node : hintPair.second)
          ss << node.toString() << " ";
        ss << endl;
      }
      return ss.str();
    }

    //serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & edges;
      ar & hints;
    }

private:
    //add edge to this ULG
    void addEdgeCore(const UnitigLinkingGraphEdge& edge, const string &seq);


    //attributes
    ISynchronizer* synchroAddEdge; //controls multithreaded access to addEdge() function
    ISynchronizer* synchroAddAllHints; //controls multithreaded access to addAllHints() function

    //store the edges of the ULG
    set<UnitigLinkingGraphEdge> edges;

    //store the hints for each trustful unitig
    map<UnitigLinkingGraphNode, set<UnitigLinkingGraphNode> > hints;
};
BOOST_CLASS_VERSION(UnitigLinkingGraph, 1)

#endif //EYTA_UNITIGLINKINGGRAPH_H
