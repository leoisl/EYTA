//
// Created by Leandro Ishi Soares de Lima on 19/04/17.
//

#ifndef KSGATB_BUBBLEOUTPUTTER_H
#define KSGATB_BUBBLEOUTPUTTER_H

#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "Utils.h"
#include "EnhanceTranscriptomeDefs.h"
#include <gatb/gatb_core.hpp>
#include <boost/algorithm/string.hpp>

#define NORMAL 0
#define DETAILED 1
#define CLUSTERED 2
#define SNP 0
#define MSNPs 1
#define REPEATS 2
#define AS 3

using namespace std;

namespace EnhanceTranscriptome {
  class BubbleOutputter {
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!
      //TODO: A BUNCH OF OUTPUT SYNCHRONIZATION IS NEEDED HERE!!!



  public:
      //constructor
      BubbleOutputter(int k, double editDistance, const string &outputPrefix, bool outputContext, const map<long, string> &transcriptIndex2sequence,
                      const map<long, set<long>> &blueUnitigId2TranscriptIndex) :
          k(k), bubbleIndex(0), editDistance(editDistance), outputContext(outputContext), synchro(System::thread().newSynchronizer()), transcriptIndex2sequence(transcriptIndex2sequence),
          blueUnitigId2TranscriptIndex(blueUnitigId2TranscriptIndex)
      {
        for (int i=NORMAL; i<=CLUSTERED;i++) {
          outputStreams.push_back(vector<ofstream*>());
          for (int j=SNP; j<=AS; j++)
            outputStreams[i].push_back(new ofstream(outputPrefix+string(".")+fileType2Description[i]+string(".")+eventType2Description[j]));
        }
      }

      ~BubbleOutputter(){
        for (int i=NORMAL; i<=CLUSTERED;i++) {
          for (int j=SNP; j<=AS; j++) {
            outputStreams[i][j]->close();
            delete outputStreams[i][j];
          }
        }
      }


      //allocate all stringstreams necessary
      void prepareOutput(const MappingInfo &mappingInfo) {
        //declare the stringstreams
        for (int i=NORMAL; i<=CLUSTERED;i++) {
          mappingInfo2StringStreams[&mappingInfo].push_back(vector<stringstream>());
          for (int j=SNP; j<=AS; j++)
            mappingInfo2StringStreams[&mappingInfo][i].push_back(stringstream());
        }
      }

      //output everything to the files and free up memory
      //TODO: check this!!
      void finalizeOutput(const MappingInfo &mappingInfo) {
        //output everything to the files NORMAL and DETAILED
        for (int i=NORMAL; i<=DETAILED;i++)
          for (int j=SNP; j<=AS; j++)
            (*outputStreams[i][j]) << mappingInfo2StringStreams[&mappingInfo][i][j].str();

        //output the clustered file
        for (int j=SNP; j<=AS; j++)
          //TODO: THIS 0.2 SHOULD BE A PARAMETER
          (*outputStreams[CLUSTERED][j]) << clusterSequences(mappingInfo2StringStreams[&mappingInfo][NORMAL][j].str(), 0.2); //sequences with 80%+ similarity will be in the same cluster


        //free up memory
        mappingInfo2AllSequencesFoundForIt.erase(&mappingInfo);
        mappingInfo2StringStreams.erase(&mappingInfo);
      }



      //checks if the new alternative path is really novel or not
      //it is novel only if it is not contained in another transcript or if it not so similar to a previous alternative path
      template <class GraphType>
      bool isAlternativePathReallyNovel (GraphType &graph, const vector<Vertex> &alternativePath, const MappingInfo &mappingInfo) {
        //build the alternative sequence
        string alternativeSequence = buildSequenceUsingULG(graph, alternativePath);

        if (checkIfAlternativeSequenceIsContainedInATranscript(graph, alternativePath, alternativeSequence)==true) {
          //this alternative path is already in a transcript - it is not novel - we should not output it
          return false;
        }

        bool alternativeSequenceIsSimilarToAPreviousOne = checkAlternativeSequenceIsSimilarToAPreviousOne(mappingInfo, alternativeSequence);

        //add this path to mappingInfo2AllSequencesFoundForIt anyway, since this is not in the transcriptome (should be output in the detailed and normal files)
        mappingInfo2AllSequencesFoundForIt[&mappingInfo].insert(alternativeSequence);
        return !alternativeSequenceIsSimilarToAPreviousOne;
      }

      //output the bubble found to a file
      //@return true if the bubble output was kind of different to other bubbles already output (at least 20% difference)
      //return false if it was similar to other bubbles already output
      template <class GraphType>
      void outputBubble(GraphType &graph, const vector<Vertex> &alternativePath,
                        const MappingInfo &mappingInfo, int startIndex, bool *firstBubble, int graphDistance, const string &prefix) {
        //get the original path without repeated nodes
        //auto posFirstNode = find(mappingInfo.getMappingFW().begin(), mappingInfo.getMappingFW().end(), alternativePath.front());
        auto posFirstNode = mappingInfo.getMappingFW().begin()+startIndex;
        auto posLastNode = find(posFirstNode+1, mappingInfo.getMappingFW().end(), alternativePath.back());
        posLastNode++;
        vector<Vertex> originalPath;
        for (auto it = posFirstNode; it!=posLastNode; ++it) {
          if (find(originalPath.begin(), originalPath.end(), *it)==originalPath.end())
            originalPath.push_back(*it);
        }


        //build the headers
        string headerOriginalSequence, headerAlternativeSequence;
        buildHeaders(headerOriginalSequence, headerAlternativeSequence, prefix);

        string headerOriginalSequenceDetailed, headerAlternativeSequenceDetailed;
        buildDetailedHeaders(headerOriginalSequenceDetailed, originalPath, mappingInfo, headerAlternativeSequenceDetailed, alternativePath, graph, prefix);

        //build the original sequence
        string originalSequence = buildSequenceUsingDBG(graph, originalPath);

        //build the alternative sequence
        string alternativeSequence = buildSequenceUsingULG(graph, alternativePath);

        //output everything
        finishOutput(mappingInfo, firstBubble, headerOriginalSequence, headerOriginalSequenceDetailed, originalSequence,
                     headerAlternativeSequence, headerAlternativeSequenceDetailed, alternativeSequence);
      }

      struct PathSorter {
          template <class GraphType>
          bool operator() (const Path<GraphType> &p1, const Path<GraphType> &p2) {
            return p1.getQuantification() > p2.getQuantification();
          }
      };


      template <class GraphType>
      void outputPathFromTheCluster(GraphType &graph, const list< Path<GraphType> > &alternativePaths, const vector<MappingInfo> &mappingInfoCluster, int cluster, bool *firstBubble) {
        /*
        //this represents the mappingInfoCluster sorted by quantification
        vector<MappingInfo> mappingInfoClusterSorted(mappingInfoCluster);
        sort(mappingInfoClusterSorted.begin(), mappingInfoClusterSorted.end(), [&](const MappingInfo &m1, const MappingInfo &m2) {
            return m1.getSRQuantification() > m2.getSRQuantification(); //from the highest to the lowest
        });

        //same thing for the alternative paths
        vector< Path<GraphType> > alternativePathsSorted(alternativePaths.begin(), alternativePaths.end());
        sort(alternativePathsSorted.begin(), alternativePathsSorted.end(), PathSorter());

        //output everything
        //get the stuff to print about the LRs into these variables
        vector<string> LROutputs(3);
        {
          stringstream LROutputSS, LROutputDetailedSS;
          bool firstLR=true;
          for (const auto &mappingInfo : mappingInfoClusterSorted) {
            string basicInfo;
            {
              stringstream ss;
              ss << ">Cluster_" << cluster
                 << (firstLR ? "_MAIN_TRANSCRIPT" : "")
                 << "_LongRead_" << mappingInfo.getReadIndex()
                 << "_Quantif_" << mappingInfo.getSRQuantification()
                 //These are just to make things easier to see in BLAT
                 << "_LR_" << mappingInfo.getReadIndex()
                 << (firstLR ? "_MT" : "");
              basicInfo = ss.str();
            }
            firstLR = false;

            //output stuff
            LROutputSS << basicInfo << endl;
            LROutputSS << mappingInfo.getRead() << endl;
            LROutputDetailedSS << basicInfo << "_Length_" << mappingInfo.getRead().size() << endl;
            LROutputDetailedSS << mappingInfo.getRead() << endl;
          }
          LROutputs[NORMAL] = LROutputSS.str();
          LROutputs[DETAILED] = LROutputDetailedSS.str();

          //TODO: put this as a parameter
          //TODO: use AI here? Unsupervised learning for clustering?
          LROutputs[CLUSTERED] = clusterSequences(LROutputs[NORMAL], 0.2); //sequences with 80%+ similarity will be in the same cluster
        }

        //declare the stringstreams
        vector<vector<stringstream> > sss;
        for (int i=NORMAL; i<=CLUSTERED;i++) {
          sss.push_back(vector<stringstream>());
          for (int j=SNP; j<=AS; j++)
            sss[i].push_back(stringstream());
        }
        vector<stringstream> onlyTheShortReads;
        for (int i=SNP; i<=AS;i++)
          onlyTheShortReads.push_back(stringstream());

        int alternativePathsIndex=0;
        for (const auto &path : alternativePathsSorted) {
          //build the alternative sequence
          string alternativeSequence = buildAlternativeSequence(graph, path.nodes);

          //find the smallest anchored distance from the alternative path to a substring of a long read
          //The idea is that the first k-1 and the last k-1 mers of the alternative paths are in the long read (since the alternative path is built from a LR node and to a LR node)
          //But we don't search for exact strings, because we want to search on all LRs
          //Thus, we get the best anchors for these k-1 mers and try to find where they fit the best
          int minAD = numeric_limits<int>::max();
          const MappingInfo* minADLR = NULL;
          for(const auto &mappingInfo : mappingInfoCluster) {
            int ad = computeAnchoredDistance(mappingInfo.getRead(), alternativeSequence, k);
            if (ad < minAD) {
              minAD = ad;
              minADLR=&mappingInfo;
            }
          }

          //get the event type
          int eventType = getEventType(minAD, alternativeSequence);

          //now, output

          //check if we have to output the LR information for this event (or if it was already output)
          if (firstBubble[eventType]) {
            for (int i=NORMAL; i<=CLUSTERED;i++)
              sss[i][eventType] << LROutputs[i];
            firstBubble[eventType]=false;
          }

          //get the basic info
          string basicInfo;
          {
            stringstream ss;
            ss << ">Cluster_" << cluster << "_AlternativePath_" << alternativePathsIndex << "_Quantif_" << path.getQuantification() << "_AP_" << alternativePathsIndex;
            basicInfo = ss.str();
          }


          //output
          sss[NORMAL][eventType] << basicInfo << endl;
          sss[NORMAL][eventType] << alternativeSequence << endl;
          sss[DETAILED][eventType] << basicInfo << "_Length_" << alternativeSequence.size()
          << "_" << eventTypeAsStr[eventType] << "_MinAnchoredDistance_" << minAD << "_with_LongRead_" << minADLR->getReadIndex() << "_Nodes=[" <<
          path.getNodesAsStr() << "]" << endl;
          sss[DETAILED][eventType] << alternativeSequence << endl;
          onlyTheShortReads[eventType] << basicInfo << endl;
          onlyTheShortReads[eventType] << alternativeSequence << endl;
          alternativePathsIndex++;
        }

        //cluster the sequences to obtain the clustered output
        for (int i=SNP; i<=AS;i++)
          //TODO: put this as a parameter?
          sss[CLUSTERED][i] << clusterSequences(onlyTheShortReads[i].str(), 0.02); //sequences with 98%+ similarity will be in the same cluster


        //output stuff concurrently
        synchro->lock();
        for (int i=NORMAL; i<=CLUSTERED;i++)
          for (int j=SNP; j<=AS; j++)
            (*outputStreams[i][j]) << sss[i][j].str();
        synchro->unlock ();
        */
      }

      //output bubbles that might contain repeats
      template <class GraphType>
      bool outputRepeatContainingBubble(GraphType &graph, const vector<Vertex> &firstPart, int numberOfNs, const vector<Vertex> &secondPart,
                                        const MappingInfo &mappingInfo, int startIndex, bool *firstBubble, int graphDistance,
                                        set< vector<Vertex> > * allOutputPaths=NULL) {
        /*
        vector<Vertex> bothParts;
        for (const auto& node : firstPart)
          bothParts.push_back(node);
        for (const auto& node : secondPart)
          bothParts.push_back(node);
        if (isItOkToOutputThisPath(bothParts, allOutputPaths)==false)
          return false;

        //build the headers
        string headerOriginalSequence, headerAlternativeSequence;
        buildHeaders(headerOriginalSequence, headerAlternativeSequence);
        string headerOriginalSequenceDetailed, headerAlternativeSequenceDetailed;
        buildDetailedHeaders(headerOriginalSequenceDetailed, headerAlternativeSequenceDetailed, mappingInfo,
                             GraphWriter::toString(firstPart, numberOfNs, secondPart, graph), graph, graphDistance);

        //build the original sequence
        string originalSequence;
        int endIndex;
        buildOriginalSequence(graph, secondPart, mappingInfo, startIndex, endIndex, originalSequence);

        //build the alternative sequence
        string alternativeSequence;
        buildAlternativeSequence(graph, firstPart, numberOfNs, secondPart, mappingInfo, startIndex, endIndex, alternativeSequence);

        //output everything
        // We lock the synchronizer
        synchro->lock();

        //output everything
        finishOutput(mappingInfo, firstBubble, headerOriginalSequence, headerOriginalSequenceDetailed, originalSequence,
                     headerAlternativeSequence, headerAlternativeSequenceDetailed, alternativeSequence);

        //increase the bubble index
        bubbleIndex++;

        // We unlock the synchronizer
        synchro->unlock ();

        return true;*/

      }

  private:
      int k;
      int bubbleIndex;
      double editDistance;
      ISynchronizer* synchro;
      static string eventTypeAsStr[];
      static string fileType2Description[];
      static string eventType2Description[];
      bool outputContext;
      vector<vector<ofstream*> > outputStreams;
      map<const MappingInfo*, set<string> > mappingInfo2AllSequencesFoundForIt;
      map<const MappingInfo*, vector< vector<stringstream> > > mappingInfo2StringStreams;
      const map<long, string> &transcriptIndex2sequence; //this will store a map to all transcripts sequences, so that we know how to search for new AS in the transcript sequences
      const map<long, set<long>> &blueUnitigId2TranscriptIndex;

      //TODO: use ML to define this?
      int getEventType(int ed, const string &alternativeSequence) const {
        double proportionOfNotMatches = (double(ed))/alternativeSequence.length();

        if (ed<=1)
          return SNP; //SNP or 1-base indel
        if (proportionOfNotMatches<=editDistance)
          return MSNPs; //Multiple SNPs / Large indels
        if (alternativeSequence.length()<=2*k-8) //TODO: change this?
          return REPEATS; //Putative repeats
        return AS; //AS
      }

      //given a string representing a fasta file:
      //>header
      //[sequence]
      //cluster the sequences and return a subset of the original sequences with copies suppressed
      //These sequences are ordered in order of importance (i.e. if the 1st, 3rd and 8th sequence falls on the same cluster, just print the 1st)
      //In my case, the order of importance is given by the quantification
      static string clusterSequences (const string &original, double threshold);
      static bool tooSimilar(const string &s, const string &t, double threshold);
      static bool tooSimilarSemiGlobal(const string &s, const string &t, double threshold);

      //check if this alternative sequence was already output or not
      bool checkAlternativeSequenceIsSimilarToAPreviousOne(const MappingInfo &mappingInfo, const string &alternativeSequence) {
          if (mappingInfo2AllSequencesFoundForIt.count(&mappingInfo)==0) {
            //no, mappingInfo2AllSequencesFoundForIt[&mappingInfo] not even exists
            return false;
          }else {
            return any_of(mappingInfo2AllSequencesFoundForIt[&mappingInfo].begin(), mappingInfo2AllSequencesFoundForIt[&mappingInfo].end(),
            [&](const string &seq) -> bool {
               //return tooSimilar(seq, alternativeSequence, 0.2);
                //TODO: THIS 0.2 SHOULD BE A PARAMETER. I AM USING THE EDIT_DISTANCE FOR NOW, BUT WE SHOULD SPLIT THESE PARAMS
                return tooSimilarSemiGlobal(seq, alternativeSequence, editDistance);
            });
          }
      }

      //check if this alternative sequence is in one of the transcripts or not
      template <class GraphType>
      bool checkIfAlternativeSequenceIsContainedInATranscript(GraphType &graph, const vector<Vertex> &alternativePath, const string &alternativeSequence) {
        //check if this path is already contained in a transcript - in this case, it is not a novel event
        //we basically check if there is a good semi-global alignment of this path to a transcript
        //the transcripts taken into account are just the ones containing any of the flanking exons
        vector<const string*> transcriptsToCheck;

        auto addTranscriptsTotranscriptsToCheck = [&](const Vertex &vertex)
        {
            long unitigId = graph[vertex].id;
            for (auto transcriptIndex : blueUnitigId2TranscriptIndex.at(unitigId))
              transcriptsToCheck.push_back(&transcriptIndex2sequence.at(transcriptIndex));
        };
        addTranscriptsTotranscriptsToCheck(alternativePath.front());
        addTranscriptsTotranscriptsToCheck(alternativePath.back());

        return any_of(transcriptsToCheck.begin(), transcriptsToCheck.end(),
                      [&](const string *seq) -> bool {
                          //TODO: THIS 0.2 SHOULD BE A PARAMETER. I AM USING THE EDIT_DISTANCE FOR NOW, BUT WE SHOULD SPLIT THESE PARAMS
                          return tooSimilarSemiGlobal(*seq, alternativeSequence, editDistance);
                      });
      }

      void buildHeaders(string& headerOriginalSequence, string& headerAlternativeSequence, const string &prefix) const {
        //build the header
        string header="";
        {
          stringstream ss;
          ss << ">Bubble_id=" << bubbleIndex << ";Length=[length];[Sequence_label]";
          header=ss.str();
        }

        //build the header of the original sequence
        headerOriginalSequence = header;
        headerOriginalSequence.replace(headerOriginalSequence.find("[Sequence_label]"),string("[Sequence_label]").length(), "Original_LR_sequence");

        //build the header of the alternative sequence
        headerAlternativeSequence = header;
        headerAlternativeSequence.replace(headerAlternativeSequence.find("[Sequence_label]"),string("[Sequence_label]").length(), prefix+"_sequence");
      }

      //build detailed headers
      template <class GraphType>
      void buildDetailedHeaders(string& headerOriginalSequence,
                                const vector<Vertex> originalPath,
                                const MappingInfo &mappingInfo,
                                string& headerAlternativeSequence,
                                const vector<Vertex> alternativePath,
                                GraphType &graph, const string &prefix) const {
        //build the header
        string header = "";
        {
          stringstream ss;
          ss << ">Bubble_id=" << bubbleIndex << ";Transcript=[transcript]" <<
          ";Length=[length];Branches=[branches];[Sequence_label];Edit_distance=[edit_distance];Read_index=" << mappingInfo.getReadIndex() <<
          ";Path=[path]";
          header = ss.str();
        }

        //add the graphDistance to the header
        /*
        string graphDistanceAsStr;
        {
          stringstream ss;
          ss << (graphDistance-1);
          graphDistanceAsStr=ss.str();
        }
        replace_first(header.find("[graph_distance]"),string("[graph_distance]").length(), graphDistanceAsStr);
         */

        //build the header of the original sequence
        headerOriginalSequence = header;
        boost::replace_first(headerOriginalSequence, "[Sequence_label]", "Original_LR_sequence");
        boost::replace_first(headerOriginalSequence, "[branches]", "N/A");
        boost::replace_first(headerOriginalSequence, "[path]", GraphWriter::toString(originalPath, graph));
        boost::replace_first(headerOriginalSequence, "[transcript]", mappingInfo.getHeader());


        //build the header of the alternative sequence
        headerAlternativeSequence = header;
        boost::replace_first(headerAlternativeSequence, "[Sequence_label]", prefix+"_sequence");
        boost::replace_first(headerAlternativeSequence, "[path]", GraphWriter::toString(alternativePath, graph));
        boost::replace_first(headerAlternativeSequence, "[transcript]", string("Novel"));

        vector<Vertex> alternativePathBranchigNodes;
        for (const auto &v : alternativePath) {
          if (isBranching(v,graph))
            alternativePathBranchigNodes.push_back(v);
        }

        boost::replace_first(headerAlternativeSequence, "[branches]", to_string(alternativePathBranchigNodes.size()));
        boost::replace_first(headerAlternativeSequence, "[branching_nodes]", GraphWriter::toString(alternativePathBranchigNodes, graph));

      }

      template <class GraphType>
      string buildSequenceUsingDBG(GraphType &graph, const vector<Vertex> &originalPath) const {
        string sequence="";
        sequence = (outputContext ? graph[originalPath[0]].seq : graph[originalPath[0]].seq.substr(graph[originalPath[0]].seq.size()-k));
        for_each(originalPath.begin()+1, originalPath.end()-1, [&](const Vertex &v) {
            sequence += graph[v].seq.substr(k - 1);
        });
        sequence += (outputContext ? graph[originalPath.back()].seq.substr(k - 1) : string(1, graph[originalPath.back()].seq[k-1]));
        return sequence;
      }


      template <class GraphType>
      string buildSequenceUsingULG(GraphType &graph, const vector<Vertex> &alternativePath) const {
        string alternativeSequence="";
        alternativeSequence = (outputContext ? graph[alternativePath[0]].seq : graph[alternativePath[0]].seq.substr(graph[alternativePath[0]].seq.size()-k));

        for(int i=1; i<alternativePath.size()-1; i++) {
          auto previousEdge = edge(alternativePath[i-1], alternativePath[i], graph).first;
          alternativeSequence += graph[previousEdge].seq;
          alternativeSequence += graph[alternativePath[i]].seq.substr(k-1);
        }

        auto previousEdge = edge(alternativePath[alternativePath.size()-2], alternativePath[alternativePath.size()-1], graph).first;
        alternativeSequence += graph[previousEdge].seq;
        alternativeSequence += (outputContext ? graph[alternativePath.back()].seq.substr(k - 1) : string(1, graph[alternativePath.back()].seq[k-1]));

        return alternativeSequence;
      }

      template <class GraphType>
      void buildAlternativeSequence(GraphType &graph, const vector<Vertex> &alternativePath,
                                    const MappingInfo &mappingInfo, int startIndex, int endIndex, string &alternativeSequence) const {
        throw runtime_error("reimplement");
        /*
        //find where to begin in the first node
        alternativeSequence = graph[alternativePath[0]].name.substr(mappingInfo.positionInsideUnitig[startIndex]);
        for_each(alternativePath.begin() + 1, alternativePath.end() - 1, [&](const Vertex &v) {
            alternativeSequence += graph[v].name.substr(k - 1);
        });
        alternativeSequence +=  graph[alternativePath.back()].name.substr(k-1,
                                                                          mappingInfo.positionInsideUnitig[endIndex]+1);
                                                                          */
      }

      template <class GraphType>
      void buildAlternativeSequence(GraphType &graph, const vector<Vertex> &firstPart, int numberOfNs, const vector<Vertex> &secondPart,
                                    const MappingInfo &mappingInfo, int startIndex, int endIndex, string &alternativeSequence) const {
        throw runtime_error("reimplement");
        /*
        //find where to begin in the first node
        alternativeSequence = graph[firstPart[0]].name.substr(mappingInfo.positionInsideUnitig[startIndex]);

        //add the whole first part
        for_each(firstPart.begin() + 1, firstPart.end(), [&](const Vertex &v) {
            alternativeSequence += graph[v].name.substr(k - 1);
        });

        //add the Ns
        for(int i=1;i<=numberOfNs;i++)
          alternativeSequence += 'N';

        //add the second part - begin in the right part
        alternativeSequence += graph[secondPart[0]].name.substr(max(0, k-1-numberOfNs));

        //add the rest of the second part
        for_each(secondPart.begin() + 1, secondPart.end() - 1, [&](const Vertex &v) {
            alternativeSequence += graph[v].name.substr(k - 1);
        });

        //add the last node
        alternativeSequence +=  graph[secondPart.back()].name.substr(k-1,
                                                                          mappingInfo.positionInsideUnitig[endIndex]+1);
                                                                          */
      }

      void finishOutput(const MappingInfo &mappingInfo, bool *firstBubble, string &headerOriginalSequence, string &headerOriginalSequenceDetailed,
                        string &originalSequence, string &headerAlternativeSequence, string &headerAlternativeSequenceDetailed, string &alternativeSequence) {
        int ed = computeEditDistance(originalSequence, alternativeSequence);

        //replace the stuff lacking in the headers
        string edAsStr;
        {
          stringstream ss;
          ss << ed;
          ss >> edAsStr;
        }
        string lengthOriginalSequenceAsStr="";
        {
          stringstream ss;
          ss << originalSequence.length();
          ss >> lengthOriginalSequenceAsStr;
        }
        string lengthAlternativeSequenceAsStr="";
        {
          stringstream ss;
          ss << alternativeSequence.length();
          ss >> lengthAlternativeSequenceAsStr;
        }
        headerOriginalSequence.replace(headerOriginalSequence.find("[length]"),string("[length]").length(),lengthOriginalSequenceAsStr);
        headerOriginalSequenceDetailed.replace(headerOriginalSequenceDetailed.find("[length]"),string("[length]").length(),lengthOriginalSequenceAsStr);
        headerOriginalSequenceDetailed.replace(headerOriginalSequenceDetailed.find("[edit_distance]"),string("[edit_distance]").length(),edAsStr);
        headerAlternativeSequence.replace(headerAlternativeSequence.find("[length]"),string("[length]").length(),lengthAlternativeSequenceAsStr);
        headerAlternativeSequenceDetailed.replace(headerAlternativeSequenceDetailed.find("[length]"),string("[length]").length(),lengthAlternativeSequenceAsStr);
        headerAlternativeSequenceDetailed.replace(headerAlternativeSequenceDetailed.find("[edit_distance]"),string("[edit_distance]").length(),edAsStr);

        //TODO: change the sequence description from Novel_AS or Novel_intron to AS or SNPs or Repeats when the file we are writing is to these
        //TODO: check the other files


        //CARE: remember to change firstBubble if the number of file increases
        int fileIndex = getEventType(ed, alternativeSequence);

        //output stuff concurrently
        synchro->lock();
        if (firstBubble[fileIndex]) {
          mappingInfo2StringStreams[&mappingInfo][NORMAL][fileIndex] << ">" << mappingInfo.getHeader() << endl;
          mappingInfo2StringStreams[&mappingInfo][NORMAL][fileIndex] << mappingInfo.getRead() << endl;
          firstBubble[fileIndex]=false;
        }

        mappingInfo2StringStreams[&mappingInfo][NORMAL][fileIndex] << headerAlternativeSequence << endl << alternativeSequence << endl;
        mappingInfo2StringStreams[&mappingInfo][DETAILED][fileIndex] << headerOriginalSequenceDetailed << endl << originalSequence << endl << headerAlternativeSequenceDetailed << endl << alternativeSequence << endl;

        //increase the bubble index
        bubbleIndex++;
        synchro->unlock ();
      }
  };
}

#endif //KSGATB_BUBBLEOUTPUTTER_H
