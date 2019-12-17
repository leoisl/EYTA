//
// Created by Leandro Ishi Soares de Lima on 18/04/17.
//

#include "EnhanceTranscriptomeDefs.h"
#include "Utils.h"

namespace EnhanceTranscriptome {
    void MappingInfo::readTheMapping(ifstream &readMappingReader, vector<Vertex> &mapping, vector<int> &positionInsideUnitig, graph_t& graph, const map<pair<int, char>, int> &labelToIndex) {
      string buffer;
      getline(readMappingReader, buffer);
      stringstream ss(buffer);
      int unitigId;
      char strand;
      int pos;
      while (ss >> unitigId >> strand >> pos) {
        if (unitigId == -1) {
          cerr << "Fatal error on MappingInfo::readTheMapping() - please report this bug" << endl;
          exit(1);
        }
        else {
          mapping.push_back(vertex(labelToIndex.at(make_pair(unitigId, strand)), graph));
          positionInsideUnitig.push_back(pos);
        }
      }
    }

    vector<MappingInfo> MappingInfo::readAllMappingInfo (const string &mappingReadsFilename, graph_t& graph, const map<pair<int, char>, int> &labelToIndex) {
      vector<MappingInfo> allMappingInfo;
      ifstream readMappingReader;
      openFileForReading(mappingReadsFilename, readMappingReader);
      string buffer;


      //read all the mapping infos
      while (getline(readMappingReader, buffer)) {
        MappingInfo mappingInfo;

        //parse the header of the mapping
        //read both indices
        {
          stringstream ss(buffer);
          int dummy;
          ss >> dummy >> mappingInfo.readIndex;
        }

        //read the header
        getline(readMappingReader, mappingInfo.header);

        //read the read
        getline(readMappingReader, mappingInfo.read);

        //read the FW mapping
        MappingInfo::readTheMapping(readMappingReader, mappingInfo.mappingFW, mappingInfo.positionInsideUnitigFW, graph, labelToIndex);

        //read the RC mapping
        MappingInfo::readTheMapping(readMappingReader, mappingInfo.mappingRC, mappingInfo.positionInsideUnitigRC, graph, labelToIndex);

        //get the unique nodes mapping here
        mappingInfo.uniqueNodes.insert(mappingInfo.mappingFW.begin(), mappingInfo.mappingFW.end());
        mappingInfo.uniqueNodes.insert(mappingInfo.mappingRC.begin(), mappingInfo.mappingRC.end());
        mappingInfo.uniqueNodes.erase(boost::graph_traits<graph_t>::null_vertex());

        //get the quantification
        mappingInfo.SRQuantification=0;
        if (mappingInfo.uniqueNodes.size()>0) {
          for (const auto &node : mappingInfo.uniqueNodes)
            mappingInfo.SRQuantification += graph[node].shortReadsCoverage;
          mappingInfo.SRQuantification /= mappingInfo.uniqueNodes.size();
        }

        //add the mapping info
        allMappingInfo.push_back(mappingInfo);
      }
      readMappingReader.close();

      return allMappingInfo;
    }

    void MappingInfo::readTranscriptSequences (const string &mappingReadsFilename, map<long, string> &transcriptIndex2sequence) {
      ifstream readMappingReader;
      openFileForReading(mappingReadsFilename, readMappingReader);
      string buffer;


      //read all the mapping infos
      while (getline(readMappingReader, buffer)) {
        //parse the header of the mapping
        long readIndex;
        string read;
        string dummy;

        //parse the header of the mapping
        //read both indices
        {
          stringstream ss(buffer);
          int dummy;
          ss >> dummy >> readIndex;
        }

        //read (and discard) the header
        getline(readMappingReader, dummy);

        //read the read
        getline(readMappingReader, read);

        //read (and discard) the FW mapping
        getline(readMappingReader, dummy);
        //read (and discard) the RC mapping
        getline(readMappingReader, dummy);

        //add it to transcriptIndex2sequence
        transcriptIndex2sequence[readIndex] = read;
      }
      readMappingReader.close();
    }
}