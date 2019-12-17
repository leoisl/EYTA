/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _TOOL_map_reads_HPP_
#define _TOOL_map_reads_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/
#include "global.h"
#include "Utils.h"
#include "UnitigCounter.h"
#include "UnitigLinkingGraph.h"
#include "build_dbg.hpp"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "MappingOutputter.h"
#include <tuple>
#include <boost/smart_ptr/scoped_ptr.hpp>
using namespace std;



////////////////////////////////////////////////////////////////////////////////
//
// THIS FILE IS AUTOMATICALLY GENERATED...
//
// THIS IS A SIMPLE EXAMPLE HOW TO USE THE Tool CLASS. IF YOU WANT MORE FEATURES,
// YOU CAN HAVE A LOOK AT THE ToyTool SNIPPET HERE:
//
//      http://gatb-core.gforge.inria.fr/snippets_tools.html
//
////////////////////////////////////////////////////////////////////////////////

//Functor class that will be cloned by the dispatcher to map the reads in a multi-threaded way
//it is a template since we might map a read or a pair of reads (in case of PE-reads mapping)
template <class SequenceOrPairOfSequence, class ReadOrFragment>
class MapAndPhase
{
protected:
    int fileIndex; //the file index of the file that we are mapping
    MappingOutputter* mappingOutputter; //class to output the results of a mapping in a synchronized way
    const string &prefix; //prefix of the files
    uint64_t &nbOfReadsProcessed; //counts the nb of reads processed
    ISynchronizer* nbOfReadsProcessedSynchro; //synchronizer for the previous var
    UnitigCounter& unitigCounter; //unitig counter
    bool printReadInTheOutput;
    const SolidKmerOracle &solidKmerOracle;
    int minSizeBlueUnitigs;
    int minSizeRedUnitigs;
    UnitigLinkingGraph* ulg;
    UnitigCounter* unitigCounterTranscripts; //required to build the ULG, to know if an unitig is blue
    map<long, set<long>>* unitigId2ReadIndex;

    //method that will be used to map the read or fragmento to the graph
    virtual void mapToTheGraph(const ReadOrFragment &header, const ReadOrFragment &readOrFragment, long readIndex) = 0;

    //method that will map a read to the graph
    void mapReadToTheGraphCore(const string &read, long readIndex, stringstream *mappingOutSS, set<UnitigInfo, UnitigInfoComparatorIdStrand> *allTrustfulUnitigsInThisRead,
                               map<long, set<long>>* unitigId2ReadIndex) {
        if (read.size() < solidKmerOracle.getKmerSize()) {
            //empty mapping
            if (mappingOutSS)
                *mappingOutSS << endl;
        }
        else {
            vector<UnitigInfo> allUnitigInfoInTheRead; allUnitigInfoInTheRead.reserve(read.size());
            set<long> allUnitigIds;

            //get all unitig info in the read
            for (int i = 0; i < read.size() - solidKmerOracle.getKmerSize() + 1; i++) {
                string LRKmer = string(read.c_str() + i, solidKmerOracle.getKmerSize());

                UnitigInfo unitigInfo;
                try {
                    //try to find the kmer in the graph (in which unitig/pos it is)
                    unitigInfo = solidKmerOracle.contains(LRKmer);

                    //kmer was found
                    stats.increaseStats("nbOfKmersInReadsThatMappedToTheGraph");
                    allUnitigIds.insert(unitigInfo.unitigId);

                    //we found the kmer in a unitig. Check if the given unitig is trustful
                    if (solidKmerOracle.isTrustableUnitig(unitigInfo.unitigId,
                                                          unitigCounterTranscripts,
                                                          minSizeBlueUnitigs, minSizeRedUnitigs)) {
                        //it is
                        stats.increaseStats("nbOfKmersInReadsThatMappedToTrustfulUnitigs");

                        //add to unitigId2ReadIndex, if it exists
                        if (unitigId2ReadIndex) {
                            (*unitigId2ReadIndex)[unitigInfo.unitigId].insert(readIndex);
                        }

                        //update allTrustfulUnitigsInThisRead, if it exists
                        if (allTrustfulUnitigsInThisRead)
                            allTrustfulUnitigsInThisRead->insert(unitigInfo);
                    }else {
                        //it is not
                        stats.increaseStats("nbOfKmersInReadsThatMappedToUntrustfulUnitigs");
                    }
                } catch (...) {
                    //kmer not found
                    unitigInfo = UnitigInfo::InvalidKmer;
                    stats.increaseStats("nbOfKmersInReadsThatDidNotMapToTheGraph");
                }
                allUnitigInfoInTheRead.push_back(unitigInfo);

                //output info for mapping
                if (mappingOutSS)
                    *mappingOutSS << unitigInfo.toString() << " ";
            }

            if (mappingOutSS)
                *mappingOutSS << endl;

            //update the unitigCounter
            unitigCounter.addToCount(allUnitigIds, 1);

            //build the ULG
            //get all pairs of trustful unitigs and build the link between them
            if (ulg) {
                vector< tuple<UnitigInfo, UnitigInfo, string> > allULGEdgesToAdd;
                vector< pair<UnitigInfo, UnitigInfo> > allHintsToAdd;
                UnitigInfo start;
                UnitigInfo end;
                string sequence;
                for (int i=0; i<allUnitigInfoInTheRead.size(); i++) {
                    //get the start unitig: the last unitig info of the next trustful unitig
                    //First, enter the next trustful unitig:
                    while (i<allUnitigInfoInTheRead.size() && solidKmerOracle.isTrustableUnitig(allUnitigInfoInTheRead[i].unitigId,
                                                                                                unitigCounterTranscripts,
                                                                                                minSizeBlueUnitigs, minSizeRedUnitigs)==false) //while I am not a trustable unitig, continue
                        i++;

                    //now, this trustful unitig is represented as a sequence of unitig info. Get the last unitig info of this trustable unitig
                    while (i+1<allUnitigInfoInTheRead.size() && allUnitigInfoInTheRead[i].unitigId == allUnitigInfoInTheRead[i+1].unitigId) //while the next unitig is equal to the previous
                        i++;

                    //start = the last unitig info of the next trustful unitig, if it exists
                    //start should be < allUnitigInfoInTheRead.size()-1 in order for it and stop to exist...
                    if (!(i<allUnitigInfoInTheRead.size()-1))
                        continue; //start or stop does not exist
                    start=allUnitigInfoInTheRead[i];

                    //pos start should always be the end of the unitig
                    //TODO: this rule-out reads with errors, where pos start could be somewhere away from the end
                    //TODO: fix this?
                    //TODO: we need a test here...
                    /*
                        if ((start.strand == 'F' && start.pos != solidKmerOracle.getUnitigSequence(start.unitigId).size()-solidKmerOracle.getKmerSize()) ||
                        (start.strand == 'R' && start.pos != 0)) {
                        */
                    if (start.pos != solidKmerOracle.getUnitigSequence(start.unitigId).size()-solidKmerOracle.getKmerSize()) {
                        //cerr << "Read = " << read << " start = " << start.toString() << " ruled out!" << endl;
                        stats.increaseStats("nbOfULGEdgesNotAddedDueToStartUnitigNotGoingUntilTheEnd");

                        //DEBUG CODE
                        /*
                        cout << "[nbOfULGEdgesNotAddedDueToStartUnitigNotGoingUntilTheEnd]: ";
                        cout << "Start: " << start.toString() << " end:" << end.toString() << " Path: ";
                        for (int k=i-5; k<i+5; k++) {
                            if (k>=0 && k<allUnitigInfoInTheRead.size())
                                cout << allUnitigInfoInTheRead[k].toString() << " ";
                        }
                        cout << endl;
                         */
                        //DEBUG CODE

                        continue; //we did not reach the end of the unitig
                    }

                    //get all end unitigs: the first unitig info of all trustful unitigs folowing the start
                    //the first end unitig will be used as a link in the ULG
                    //all end unitigs are used to provide hints
                    bool firstEndUnitig=true;
                    for (int j=i+1; //start on the next unitig info
                         j<allUnitigInfoInTheRead.size();
                         j++) {

                        //go to the next trustful unitig
                        while (j<allUnitigInfoInTheRead.size() && solidKmerOracle.isTrustableUnitig(allUnitigInfoInTheRead[j].unitigId,
                                                                                                    unitigCounterTranscripts,
                                                                                                    minSizeBlueUnitigs, minSizeRedUnitigs)==false) //while I am not a trustable unitig, continue
                            j++;

                        //check if end exists
                        if (!(j<allUnitigInfoInTheRead.size()))
                            continue; //end does not exist

                        //I entered here the first position of the next trustful unitig!
                        end=allUnitigInfoInTheRead[j];

                        if (start.unitigId == end.unitigId) {
                            //if we are linking 2 equal unitigs due to missing kmers in the short reads , skip
                            stats.increaseStats("nbOfULGEdgesNotAddedDueToMappingToEqualUnitigs");

                            //but before, go until the end of the unitig
                            //we go to the last unitig info of this trustful unitig, in order to get the next one in the next loop
                            while (j+1<allUnitigInfoInTheRead.size() && allUnitigInfoInTheRead[j].unitigId == allUnitigInfoInTheRead[j+1].unitigId) //while the next unitig is equal to the previous
                                j++;

                            continue;
                        }

                        //pos end should always be the begin of the unitig
                        //TODO: this rule-out reads with errors, where pos end could be somewhere away from the begin
                        //TODO: fix this?
                        //TODO: test this
                        /*
                        if ((end.strand == 'F' && end.pos != 0) ||
                            (end.strand == 'R' && end.pos != solidKmerOracle.getUnitigSequence(end.unitigId).size()-solidKmerOracle.getKmerSize())) {
                        */
                        if (end.pos != 0) {
                            //cerr << "Read = " << read << " end = " << end.toString() << " ruled out!" << endl;

                            if (firstEndUnitig) {
                                stats.increaseStats("nbOfULGEdgesNotAddedDueToEndUnitigNotStartingInTheBegin");

                                //DEBUG CODE
                                /*
                                cout << "[nbOfULGEdgesNotAddedDueToEndUnitigNotStartingInTheBegin]: ";
                                cout << "Start: " << start.toString() << " end:" << end.toString() << " Path: ";
                                for (int k=i-5; k<j+5; k++) {
                                    if (k>=0 && k<allUnitigInfoInTheRead.size())
                                        cout << allUnitigInfoInTheRead[k].toString() << " ";
                                }
                                cout << endl;
                                 */
                                //DEBUG CODE
                            }
                            else {
                                stats.increaseStats("nbOfHintsNotAddedDueToEndUnitigNotStartingInTheBegin");

                                //DEBUG CODE
                                /*
                                cout << "[nbOfHintsNotAddedDueToEndUnitigNotStartingInTheBegin]: ";
                                cout << "Start: " << start.toString() << " end:" << end.toString() << " Path: ";
                                for (int k=i-5; k<j+5; k++) {
                                    if (k>=0 && k<allUnitigInfoInTheRead.size())
                                        cout << allUnitigInfoInTheRead[k].toString() << " ";
                                }
                                cout << endl;
                                 */
                                //DEBUG CODE
                            }


                            //but before, go until the end of the unitig
                            //we go to the last unitig info of this trustful unitig, in order to get the next one in the next loop
                            while (j+1<allUnitigInfoInTheRead.size() && allUnitigInfoInTheRead[j].unitigId == allUnitigInfoInTheRead[j+1].unitigId) //while the next unitig is equal to the previous
                                j++;

                            continue; //we did not start at the begin of the unitig
                        }

                        //create the ULG edge if the start and end are consecutive trustful unitig
                        if (firstEndUnitig) {
                            firstEndUnitig=false;

                            //I get the sequence between the start and the end to be the label of the arc between them
                            sequence = read.substr(i+solidKmerOracle.getKmerSize(), j-i-1);

                            //We save the edge to add it
                            allULGEdgesToAdd.push_back(make_tuple(start, end, sequence));

                            //add it
                            stats.increaseStats("nbOfULGEdgesAdded");
                        }

                        //add this unitig as a hint anyway, since start can reach it
                        allHintsToAdd.push_back(make_pair(start, end));
                        stats.increaseStats("nbOfHintsAdded");


                        //we go to the last unitig info of this trustful unitig, in order to get the next one in the next loop
                        while (j+1<allUnitigInfoInTheRead.size() && allUnitigInfoInTheRead[j].unitigId == allUnitigInfoInTheRead[j+1].unitigId) //while the next unitig is equal to the previous
                            j++;
                    }
                }

                //add all the edges saved in allULGEdgesToAdd
                ulg->addEdges(allULGEdgesToAdd, solidKmerOracle);

                //add all hints found
                ulg->addAllHints(allHintsToAdd, solidKmerOracle);
            }
        }
    }

public:
    MapAndPhase (int fileIndex, MappingOutputter* mappingOutputter, const string &prefix, uint64_t &nbOfReadsProcessed, ISynchronizer* nbOfReadsProcessedSynchro, UnitigCounter &unitigCounter,
                 bool printReadInTheOutput, const SolidKmerOracle &solidKmerOracle,
                 int minSizeBlueUnitigs, int minSizeRedUnitigs, UnitigLinkingGraph* ulg, UnitigCounter* unitigCounterTranscripts, map<long, set<long>>* unitigId2ReadIndex) :
        fileIndex(fileIndex), mappingOutputter(mappingOutputter), prefix(prefix),
        nbOfReadsProcessed(nbOfReadsProcessed), nbOfReadsProcessedSynchro(nbOfReadsProcessedSynchro), unitigCounter(unitigCounter),
        printReadInTheOutput(printReadInTheOutput),solidKmerOracle(solidKmerOracle),
        minSizeBlueUnitigs(minSizeBlueUnitigs), minSizeRedUnitigs(minSizeRedUnitigs), ulg(ulg), unitigCounterTranscripts(unitigCounterTranscripts), unitigId2ReadIndex(unitigId2ReadIndex){}

    //method that will be called for each read / pair of read
    virtual void operator()(const SequenceOrPairOfSequence& sequenceOrPairOfSequence) = 0;
};

//specialization of MapAndPhase for single-end reads
class MapAndPhaseSERead : public MapAndPhase<Sequence, string> {
protected:
    virtual void mapToTheGraph(const string& header, const string &read, long readIndex) {
        //this will store the mapping output of this read
        stringstream mappingOutSS;

        //mapping
        if (mappingOutputter){
            mappingOutSS << fileIndex << " " << readIndex;
            if (printReadInTheOutput)
                mappingOutSS << endl << header << endl << read;
            mappingOutSS << endl;
        }

        //map the read
        mapReadToTheGraphCore(read, readIndex, &mappingOutSS, NULL, unitigId2ReadIndex);

        //map the RC
        //we need to map this, unless the data is stranded
        mapReadToTheGraphCore(reverse_complement(read), readIndex, &mappingOutSS, NULL, unitigId2ReadIndex);

        //output the mapping found
        if (mappingOutputter)
            mappingOutputter->output(mappingOutSS.str());
    }

public:
    using MapAndPhase<Sequence, string>::MapAndPhase; //inherit the constructor from base class

    //template-specific method overloading
    virtual void operator()(const Sequence& seq) {
        string read = seq.toString();
        read = toUpperCase(read);

        //map this read to the graph
        mapToTheGraph(seq.getComment(), read, seq.getIndex());
    }
};

//specialization of MapAndPhase for paired-end reads
class MapAndPhasePERead : public MapAndPhase< pair<Sequence, Sequence>, pair<string, string> > {
protected:
    virtual void mapToTheGraph(const pair<string, string> &pairOfHeaders, const pair<string, string> &pairOfReads, long readIndex) {
        //this will store the mapping output
        stringstream mappingOutSS;
        string read;

        //map the fragment in a FW-way
        {
            //map the FW left read
            read = pairOfReads.first;
            if (mappingOutputter) {
                mappingOutSS << fileIndex << " " << readIndex;
                if (printReadInTheOutput)
                    mappingOutSS << " " << read;
                mappingOutSS << endl;
            }
            set<UnitigInfo, UnitigInfoComparatorIdStrand> allTrustfulUnitigsInFWLeftRead;
            mapReadToTheGraphCore(read, readIndex, &mappingOutSS, &allTrustfulUnitigsInFWLeftRead, NULL);

            //map the RC right read
            read = reverse_complement(pairOfReads.second);
            if (mappingOutputter){
                mappingOutSS << (fileIndex+1) << " " << readIndex;
                if (printReadInTheOutput)
                    mappingOutSS << " " << read;
                mappingOutSS << endl;
            }
            set<UnitigInfo, UnitigInfoComparatorIdStrand> allTrustfulUnitigsInRCRightRead;
            mapReadToTheGraphCore(read, readIndex, &mappingOutSS, &allTrustfulUnitigsInRCRightRead, NULL);

            //create the hints from allTrustfulUnitigsInFWLeftRead to allTrustfulUnitigsInRCRightRead
            {
                vector<pair<UnitigInfo, UnitigInfo> > allHintsToAdd;
                for (const UnitigInfo &leftUnitig : allTrustfulUnitigsInFWLeftRead) {
                    for (const UnitigInfo &rightUnitig : allTrustfulUnitigsInRCRightRead) {
                        allHintsToAdd.push_back(make_pair(leftUnitig, rightUnitig));
                    }
                }
                //add all hints
                ulg->addAllHints(allHintsToAdd, solidKmerOracle);
            }
        }

        //map the fragment in a RC-way
        {
            //map the FW right read
            read = pairOfReads.second;
            if (mappingOutputter){
                mappingOutSS << (fileIndex+1) << " " << readIndex;
                if (printReadInTheOutput)
                    mappingOutSS << " " << read;
                mappingOutSS << endl;
            }
            set<UnitigInfo, UnitigInfoComparatorIdStrand> allTrustfulUnitigsInFWRightRead;
            mapReadToTheGraphCore(read, readIndex, &mappingOutSS, &allTrustfulUnitigsInFWRightRead, NULL);

            //map the RC left read
            read = reverse_complement(pairOfReads.first);
            if (mappingOutputter){
                mappingOutSS << fileIndex << " " << readIndex;
                if (printReadInTheOutput)
                    mappingOutSS << " " << read;
                mappingOutSS << endl;
            }
            set<UnitigInfo, UnitigInfoComparatorIdStrand> allTrustfulUnitigsInRCLeftRead;
            mapReadToTheGraphCore(read, readIndex, &mappingOutSS, &allTrustfulUnitigsInRCLeftRead, NULL);

            //create the hints from allTrustfulUnitigsInFWRightRead to allTrustfulUnitigsInRCLeftRead
            {
                vector<pair<UnitigInfo, UnitigInfo> > allHintsToAdd;
                for (const UnitigInfo &leftUnitig : allTrustfulUnitigsInFWRightRead) {
                    for (const UnitigInfo &rightUnitig : allTrustfulUnitigsInRCLeftRead) {
                        allHintsToAdd.push_back(make_pair(leftUnitig, rightUnitig));
                    }
                }
                //add all hints
                ulg->addAllHints(allHintsToAdd, solidKmerOracle);
            }
        }


        //output the mapping found
        if (mappingOutputter)
            mappingOutputter->output(mappingOutSS.str());
    }

public:
    using MapAndPhase<pair<Sequence, Sequence>, pair<string, string>>::MapAndPhase; //inherit the constructor from base class

    //template-specific method overloading
    virtual void operator()(const pair<Sequence, Sequence> & pairSeq) {
        string leftRead = pairSeq.first.toString();
        leftRead = toUpperCase(leftRead);
        string rightRead = pairSeq.second.toString();
        rightRead = toUpperCase(rightRead);

        //map this pair of reads to the graph
        mapToTheGraph(make_pair(pairSeq.first.getComment(), pairSeq.second.getComment()),
                      make_pair(leftRead, rightRead), pairSeq.first.getIndex());
    }
};

class map_reads : public Tool
{
private:
    void checkParameters();

public:

    // Constructor
    map_reads ();

    // Actual job done by the tool is here
    void execute ();
};

/********************************************************************************/

#endif
