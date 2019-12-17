//
// Created by Leandro Ishi Soares de Lima on 10/10/17.
//

#include "debug.h"

#ifdef EYTA_DEBUG
//DEBUG
void printKmersToFile (const Graph& graph, string filePath, string edgesPath) {
  ofstream file(filePath);
  file << "node\tseq\tabundance" << endl;
  ofstream fileEdges(edgesPath);
  fileEdges << "source\ttarget" << endl;

  // We get an iterator for all nodes of the graph.
  GraphIterator<Node> it = graph.iterator();

  // We loop each node. Note the structure of the for loop.
  for (it.first(); !it.isDone(); it.next())
  {
    // We dump an ascii representation of the current node.
    file << graph.nodeMPHFIndex (it.item()) << "\t" << graph.toString (it.item()) << "\t" << graph.queryAbundance(it.item()) << std::endl;


    graph.successors(it.item()).iterate([&](Node &node) {
        fileEdges << graph.nodeMPHFIndex (it.item()) << "\t" << graph.nodeMPHFIndex (node) << endl;
    });

    graph.predecessors(it.item()).iterate([&](Node &node) {
        fileEdges << graph.nodeMPHFIndex (node) << "\t" << graph.nodeMPHFIndex (it.item()) << endl;
    });
  }
  file.close();
  fileEdges.close();
}


//create the solid_kmers_binary_with_count file: a file with 2 columns, where the first is a k-mer encoded in 2-bit format and the second is the kmer abundance
void dump_solid_kmers(const gatb::core::debruijn::impl::Graph& graph, const string& kmer_counts_name)
{
  //create the output file
  ofstream kmer_counts_file(kmer_counts_name, ofstream::binary);

  //create an iterator of the nodes
  ProgressGraphIterator<Node, ProgressTimer> it (
      graph.iterator(), "Graph: dumping solid kmers with count");

  //output all nodes
  bool firstKmer=true;
  int kmerSize;
  int kmerNBits;
  for (it.first(); !it.isDone(); it.next()) {
    if (firstKmer) {
      kmerSize = (int)(graph.getKmerSize());
      kmerNBits = (int)(it.item().kmer.getSize());

      // write how many bits the type storing the kmer has (8, 16, 24, etc...)
      kmer_counts_file.write((const char*)&kmerNBits, sizeof(kmerNBits));

      //write k-mer size as the next 4 bits (an int)
      kmer_counts_file.write((const char*)&kmerSize, sizeof(kmerNBits));

      firstKmer = false;
    }

    //write a kmer
    gatb::core::tools::math::Integer &kmer = it.item().kmer;
    int writtenBytes = 0;

    for (int i = 0; i < kmerSize; ) {
      u_int8_t byte = 0;

      for (int j = 0; j < 4; j++, i++) {
        if (i>=kmerSize)
          break;
        byte |= (kmer[i] << (j * 2));
      }
      kmer_counts_file.write((const char *) &byte, sizeof(u_int8_t));
      ++writtenBytes;
    }

    //finish the writing
    u_int8_t zeroByte=0;
    for (;writtenBytes<kmerNBits/8;++writtenBytes)
      kmer_counts_file.write((const char *) &zeroByte, sizeof(u_int8_t));

    //write the abundance
    int abundanceAsInt = (int) (it.item().abundance);
    kmer_counts_file.write((const char *) &abundanceAsInt, sizeof(int));
  }

  kmer_counts_file.close();
}
#endif