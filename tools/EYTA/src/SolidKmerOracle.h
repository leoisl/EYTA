//
// Created by Leandro Ishi Soares de Lima on 10/10/17.
//

#ifndef KSGATB_SOLIDKMERORACLE_H
#define KSGATB_SOLIDKMERORACLE_H

#include <vector>
#include <string>
#include <sstream>
#include <gatb/gatb_core.hpp>
#include "Tests.h"
#include "ExceptionWithFileAndLine.h"
#include "UnitigCounter.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

using namespace std;


//represents the unitig id a GATB node maps to
//pos is not always defined - it is optional
struct UnitigInfo {
    long unitigId;
    char strand;
    int pos;
    static UnitigInfo InvalidKmer; //represents an invalid kmer

    UnitigInfo (long unitigId, char strand, int pos=-1) : unitigId(unitigId), strand(strand), pos(pos){}
    UnitigInfo(){}
    string toString() const;
    bool isValid() {return unitigId>=0; }
};

//unitig info comparator with id and strand - used in map_reads.hpp
struct UnitigInfoComparatorIdStrand {
    bool operator()(const UnitigInfo &u1, const UnitigInfo& u2) const {
        if (u1.unitigId != u2.unitigId)
            return u1.unitigId < u2.unitigId;
        return u1.strand < u2.strand;
    }
};

//exceptions in case of errors
class InexistentKmer : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};
class InvalidUnitig : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};
class InvalidKmer : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};

//Used to store the unitigs sequence and to check if a kmer really exist in the graph
class SolidKmerOracle {
public:
    //constructor
    //Graph: the graph to work on
    SolidKmerOracle(Graph* graph, const string &graphPrefix) : graph(graph), graphPrefix(graphPrefix)
        ,unitigsSequences(), nodeIdToUnitigId(((size_t)graph->getInfo()["kmers_nb_solid"]->getInt())) {}

    //constructor strictly for the serialization
    SolidKmerOracle(){}
    ~SolidKmerOracle(){
        delete graph;
    }

    //@brief: Add an unitig sequence that is present on the graph, by adding the sequence to unitigsSequences and also filling up nodeIdToUnitigId
    //@params:
    //  seq: the unitig sequence - this is intentionally seq and not const string &seq
    //  startingNode: the first node of the unitig sequence in order to traverse the graph
    void addUnitigSequence (string seq, const Node &leftestNode);

    //@brief: checks if the kmer given as input is in the graph or not. If it is, return an UnitigIdStrandPos of the kmer. Otherwise, throws either an invalid kmer or an inexistent kmer exception
    //this is intentionally string kmer and not const string &kmer
    UnitigInfo contains (string kmer) const;

    //get the unitig sequence stored in this oracle
    string getUnitigSequence (long index) const{
        if (index>=unitigsSequences.size())
            throw InvalidUnitig(__FILE__, __LINE__);
        return unitigsSequences[index];
    }

    bool isTrustableUnitig (long index, UnitigCounter* unitigCounterTranscripts, int minSizeBlueUnitigs, int minSizeRedUnitigs) const {
        if (index==-1) return false; //invalid unitigs are never trustable
        if (index>=unitigsSequences.size()) return false; //unitigs that are not in the range are never trustable
        int size = getUnitigSequence(index).size(); //get the size of the unitig

        if (unitigCounterTranscripts==NULL)
            //here, the unitigCounterTranscripts is NULL... we suppose then we are working with blue unitigs
            return size>=minSizeBlueUnitigs;


        if (unitigCounterTranscripts->getNodeCount(index)>0 && size>=minSizeBlueUnitigs) return true; //if the unitig is blue, check if it has the minimum size
        return size >= minSizeRedUnitigs; //check if it satifies the minimim red size
    }

    //get the nb of unitigs
    long getNbOfUnitigs () const {
        return unitigsSequences.size();
    }

    //get the kmer size
    long getKmerSize() const {
        return graph->getKmerSize();
    }
private:
    //for easy testing
    friend class EYTATester;

    //attributes
    Graph* graph; //the graph where this oracle works on
    string graphPrefix; //graph is written to graphPrefix+".h5"

    //TODO: vector<string> is not good here... when it reaches the limit, it is going to double and copy all the elements...
    //TODO: find something more efficient...
    vector<string> unitigsSequences; //the unitig sequences themselves - TODO: we should not use string, but a bit array (look at codeseed of the ModelDirect, Canonical, etc)
    //TODO: find something more efficient...

    vector < long > nodeIdToUnitigId; //maps node to unitig id, strand pos

    //serialization
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar  & graphPrefix;
        ar  & unitigsSequences;
        ar  & nodeIdToUnitigId;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar  & graphPrefix;
        ar  & unitigsSequences;
        ar  & nodeIdToUnitigId;

        //when loading the object, load also the graph
        graph = gatb::core::debruijn::impl::Graph::loadAsPointer(graphPrefix+".h5");
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};
BOOST_CLASS_VERSION(SolidKmerOracle, 1)

#endif //KSGATB_SOLIDKMERORACLE_H
