//
// Created by Leandro Ishi Soares de Lima on 10/10/17.
//

#include "SolidKmerOracle.h"
#include "Utils.h"
#include <boost/algorithm/string/predicate.hpp>

UnitigInfo UnitigInfo::InvalidKmer(-1, '?', -1);

string UnitigInfo::toString() const {
  stringstream ss;
  ss << unitigId << " " << strand << " " << pos;
  return ss.str();
}

//@brief: Add an unitig sequence that is present on the graph, by adding the sequence to unitigsSequences and also filling up nodeIdToUnitigId
//@params:
//  seq: the unitig sequence - this is intentionally seq and not const string &seq
//  startingNode: the first node of the unitig sequence in order to traverse the graph
void SolidKmerOracle::addUnitigSequence (string seq, const Node &leftestNode) {
  seq = toUpperCase(seq);

  //maps a char to a Nucleotide
  map<char,Nucleotide> char2nt;
  char2nt['A']=NUCL_A;
  char2nt['C']=NUCL_C;
  char2nt['T']=NUCL_T;
  char2nt['G']=NUCL_G;

  //populate nodeIdToUnitigId with the new unitig id given
  //associate the node id to its unitig id
  //Note: GATB kmer is any kmer... It is not the canonical one (i.e. smaller one). Maybe is the one that was added...
  //Anyway, for a given kmer, GATB define one of them as the forward node. This forward node is what it counts (and it is not necessarily the canonical kmer)
  auto currentNode = leftestNode;
  int i=0;

  for_each(seq.begin() + graph->getKmerSize() - 1, seq.end(), [&](char c) {
      if (char2nt.find(c)==char2nt.end()) {
        cout << "Fatal bug on SolidKmerOracle::addUnitigSequence(): invalid character: " << c << endl;
        cout.flush();
        exit(1);
      }

      Nucleotide nucleotide = char2nt[c];

      if (i>0)
        currentNode = graph->successor(currentNode, nucleotide);

      nodeIdToUnitigId[graph->nodeMPHFIndex(currentNode)] = unitigsSequences.size();
      i++;
  });

  //add the unitig sequence itself
  unitigsSequences.push_back(seq);
}

//@brief: checks if the kmer given as input is in the graph or not. If it is, return an UnitigIdStrandPos of the kmer. Otherwise, throws either an invalid kmer or an inexistent kmer exception
//this is intentionally string kmer and not const string &kmer
UnitigInfo SolidKmerOracle::contains (string kmer) const {
  kmer = toUpperCase(kmer);

  //check the kmer size
  if (kmer.size() != graph->getKmerSize())
    throw InvalidKmer(__FILE__, __LINE__);

  //TODO: From my tests, if you use buildNode with a Kmer containing an 'N', it will build a node with all Ns replaced by G
  //TODO: However, Kmers containing Ns are NOT included in any unitig (i.e. the graph builder does not convert an N to a G, and build the graph. It simply disregards kmers containing Ns)
  //TODO: this is a sanity check to also discard Kmers containing Ns
  //TODO: check this with GATB team

  //TODO: UPDATE
  //TODO: Magali had a dataset where we had the base 'K' in the fasta file
  //TODO: So I am just discarding all reads that are not composed by ACGT
  if (!boost::all(kmer, [](char c) -> bool {
      return c=='A' || c=='a' || c=='C' || c=='c' || c=='G' || c=='g' ||c=='T' || c=='t';
  })) {
      throw InvalidKmer(__FILE__, __LINE__);
  }

  //build the node
  Node node = graph->buildNode(kmer.c_str());

  //get the unitig localization of this kmer
  u_int64_t index = graph->nodeMPHFIndex(node);

  //check if the index is valid
  if (index>=nodeIdToUnitigId.size())
    throw InexistentKmer(__FILE__, __LINE__);

  //get the unitig info - this is where the GATB node maps to the unitigs
  long unitigId=nodeIdToUnitigId[index];

  //get the strand and position
  char strand;
  size_t pos;
  string unitigSequence = getUnitigSequence(unitigId);
  if ((pos=unitigSequence.find(kmer))!=string::npos) {
    //the kmer and the sequence of the unitig match
    strand='F';
  }else if ((pos=reverse_complement(unitigSequence).find(kmer))!=string::npos) {
    strand='R';
  }else {
    //the kmer does not match the unitig: it is a false positive
    //kmer not found
    throw InexistentKmer(__FILE__, __LINE__);
  }

  //kmer found, return the result
  return UnitigInfo(unitigId, strand, (int)pos);
}

