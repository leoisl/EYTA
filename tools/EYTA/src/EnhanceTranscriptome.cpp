//
// Created by Leandro Ishi Soares de Lima on 27/01/17.
//

#include "EnhanceTranscriptome.h"
#include "global.h"
#include "Utils.h"
#include <gatb/gatb_core.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace EnhanceTranscriptome {
    //constructor
    EnhanceTranscriptome::EnhanceTranscriptome(const string &prefix,
                                               const string &outputPrefix, int k,
                                               int maxLengthAS, int maxLengthIntron, double editDistance,
                                               int splicingComplexityAS,
                                               int splicingComplexityIR,
                                               bool outputContext,
                                               int nbCores) :
          graph(countNbNodes(prefix+".unitigs.size")*2+1), //*2 because of the RC, and +1 for the artificial node
          k(k),
          maxLengthAS(maxLengthAS+1), //1 because there is an off by 1
          maxLengthIntron(maxLengthIntron+1), //1 because there is an off by 1
          splicingComplexityAS(splicingComplexityAS),
          splicingComplexityIR(splicingComplexityIR)
    {
      //create the nodes of the graph
      string nodesFilename(prefix+".colored.nodes");
      cerr << "Creating the nodes of the Boost graph..." << endl;
      {
        ifstream nodesFileReader(nodesFilename);
        int id;
        int isTrustable;
        int transcriptomicCount;
        int shortReadsCount;
        string seq;
        int index = 0;
        while (nodesFileReader >> id >> isTrustable >> transcriptomicCount >> shortReadsCount >> seq) {
          //add the FW node
          Vertex vF = vertex(index, graph);
          labelToIndex[make_pair(id, 'F')] = index;
          indexToLabel[index] = make_pair(id, 'F');

          //TODO: sequence should be stored as a 2-bit array
          graph[vF].id = id;
          graph[vF].isTrustable = bool(isTrustable);
          graph[vF].transcriptomeCoverage = transcriptomicCount;
          graph[vF].shortReadsCoverage = shortReadsCount;
          graph[vF].seq = seq;
          graph[vF].strand = 'F';
          graph[vF].weight = seq.length()-k+1;
          index++;

          //add the RC node
          Vertex vR = vertex(index, graph);
          labelToIndex[make_pair(id, 'R')] = index;
          indexToLabel[index] = make_pair(id, 'R');

          //TODO: sequence should be stored as a 2-bit array
          graph[vR].id = id;
          graph[vR].isTrustable = bool(isTrustable);
          graph[vR].transcriptomeCoverage = transcriptomicCount;
          graph[vR].shortReadsCoverage = shortReadsCount;
          graph[vR].seq = reverse_complement(seq);
          graph[vR].strand = 'R';
          graph[vR].weight = seq.length()-k+1;
          index++;
        }
        nodesFileReader.close();

        //build the artificial node
        //TODO: do we need this?
        artificialIndex=index;
        Vertex artificial = vertex(artificialIndex, graph);
        labelToIndex[make_pair(ARTIFICIAL_ID, 'F')] = artificialIndex;
        indexToLabel[artificialIndex] = make_pair(ARTIFICIAL_ID, 'F');
        graph[artificial].id = ARTIFICIAL_ID;
        graph[artificial].isTrustable = true;
        graph[artificial].transcriptomeCoverage = 0;
        graph[artificial].shortReadsCoverage = 0;
        graph[artificial].seq = "";
        graph[artificial].strand = 'F';
        graph[artificial].weight = 0;
      }
      cerr << "Creating the nodes of the Boost graph... - Done!" << endl;




      //add the edges from the ULG
      string ulgFilename(prefix+".ulg");
      //create the edges of the graph - these edges are based on the ULG edges
      cerr << "Creating the edges of the Boost graph based on the " << ulgFilename << "..." << endl;
      UnitigLinkingGraph ulg;
      cerr << "Loading " << ulgFilename << " ..." << endl;
      std::ifstream ulgIF(ulgFilename);
      boost::archive::text_iarchive textIArchive(ulgIF);
      textIArchive >> ulg;
      cerr << "Loading " << ulgFilename << " ... - Done!" << endl;
      {
        int index = 0;
        for (const auto &edge : ulg.getEdges()) {
          //tries to add the edge
          pair<Edge, bool> return_from_add_edge = add_edge(vertex(labelToIndex[make_pair(edge.getFrom().getId(), edge.getFrom().getStrand())], graph),
                                                           vertex(labelToIndex[make_pair(edge.getTo().getId(), edge.getTo().getStrand())], graph),
                                                           graph);
          if (return_from_add_edge.second) {
            //if it was sucessfully added, configure it
            graph[return_from_add_edge.first].id = index;
            //TODO: IF WE ARE GOING TO USE THE SEQUENCE CONTENT, WE SHOULD USE ALL SEQUENCES, NOT ONLY THE REPRESENTATIVE...
            graph[return_from_add_edge.first].seq = edge.getRepresentativeSequence();
            graph[return_from_add_edge.first].edgeWeight = edge.getRepresentativeSequence().size();

            //the target weight is the weight of the target node + the weight of this edge, because to reach the target, you have to go through this edge anyway
            graph[return_from_add_edge.first].targetWeight = graph[target(return_from_add_edge.first, graph)].weight + graph[return_from_add_edge.first].edgeWeight;
            graph[return_from_add_edge.first].sourceWeight = graph[source(return_from_add_edge.first, graph)].weight + graph[return_from_add_edge.first].edgeWeight;

            index++;
          }
        }
      }
      cerr << "Creating the edges of the Boost graph... - Done!" << endl;


      //EYTA_DEBUG
      //print the graph to be viewed in cytoscape
      //GraphWriter::writeGraphToFile(graph, "graph.edges.cytoscape");
      //EYTA_DEBUG


      //load the whole transcriptome in memory
      cerr << "Loading transcriptome..." << endl;
      string mappingReadsFilename(prefix+".tr.mapping.0");
      map<long, string> transcriptIndex2sequence;
      MappingInfo::readTranscriptSequences(mappingReadsFilename, transcriptIndex2sequence);
      cerr << "Loading transcriptome... - Done!" << endl;

      //load blue_trustable_nodes_to_transcripts
      cerr << "Loading " << (prefix+"_blue_trustable_nodes_to_transcripts")  << "..." << endl;
      //this will associate each blue unitig to the transcripts it is contained in
      map<long, set<long>> blueUnitigId2TranscriptIndex;

      {
        ifstream blueTrustableNodesToTranscriptsFile(prefix + "_blue_trustable_nodes_to_transcripts");
        boost::archive::text_iarchive textIArchive(blueTrustableNodesToTranscriptsFile);
        textIArchive >> blueUnitigId2TranscriptIndex;
      }
      cerr << "Loading " << (prefix+"_blue_trustable_nodes_to_transcripts")  << "... - Done!" << endl;



      cerr << "[Starting Enhancing Transcriptome process... ]" << endl;

      {
        //declare the graph filters and the filtered graph
        SubgraphEdgeFilter* subgraphEdgeFilter = new SubgraphEdgeFilter(&graph);
        SubgraphVertexFilter* subgraphVertexFilter = new SubgraphVertexFilter(&graph, subgraphEdgeFilter);
        FilteredGraph* filteredGraph = new FilteredGraph(graph, SubgraphEdgeFilterForFG<SubgraphEdgeFilter>(subgraphEdgeFilter),
                                          SubgraphVertexFilterForFG<SubgraphVertexFilter>(subgraphVertexFilter));


        //create the bubble outputter
        BubbleOutputter bubbleOutputter(k, editDistance, outputPrefix, outputContext, transcriptIndex2sequence, blueUnitigId2TranscriptIndex);

        //read the mapping
        vector<MappingInfo> allMappingInfos = MappingInfo::readAllMappingInfo(mappingReadsFilename, graph, labelToIndex);


        //TODO: change this 0
        auto strategyImpl = Strategy::getStrategy(0, num_vertices(graph), ulg, labelToIndex);

        //process all transcripts
        long nbOfTranscriptsProcessed = 0;
        long timeSpendProcessingInSeconds=0;

        //TODO: MULTITHREADED PART MUST BE IMPLEMENTED HERE!!!
        for (const auto &mappingInfo : allMappingInfos) {
          cerr << "Processing transcript " << mappingInfo.getReadIndex() << endl;
          //time measuring
          boost::posix_time::ptime startTime(boost::posix_time::microsec_clock::local_time());

          //put everyone in the graph
          subgraphVertexFilter->reset(true);

          //remove all blue-only nodes - we don't need 'em - we are searching for red paths
          //and we do not want to risk to assembling chimeric paths by going to regions that belong to other genes
          //nor searching for events that are described by other transcripts in the same gene
          subgraphVertexFilter->removeAllBlueOnlyNodes();

          //process the transcript
          strategyImpl->processTranscript(mappingInfo, *filteredGraph, graph, subgraphVertexFilter, maxLengthAS,
                                          maxLengthIntron,
                                          splicingComplexityAS, splicingComplexityIR, k, bubbleOutputter, ulg, k - 10,
                                          k + 1);


          cout << "Finished processing transcript " << mappingInfo.getReadIndex() << endl;
          boost::posix_time::ptime endTime(boost::posix_time::microsec_clock::local_time());
          long secondsTaken = (endTime - startTime).total_seconds();
          timeSpendProcessingInSeconds += secondsTaken;
          nbOfTranscriptsProcessed++;
          double averageProcessingTime = ((double) (timeSpendProcessingInSeconds)) / nbOfTranscriptsProcessed;
          cout << "Average seconds for processing a transcript: " << (long) averageProcessingTime << endl;
          //TODO: Removed this
          //cout << "Estimated remaining time (in seconds): " << ((long)((nbOfMappingOutputFiles-nbOfTranscriptsProcessed)*averageProcessingTime/nbCores)) << endl;
        }
      }

      cerr << endl << "[Enhancing Transcriptome process finished!]" << endl;
      cerr.flush();
    }
}