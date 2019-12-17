//
// Created by Leandro Ishi Soares de Lima on 10/10/17.
//

#ifndef KSGATB_DEBUG_H
#define KSGATB_DEBUG_H
#include "global.h"
#include <gatb/gatb_core.hpp>

#ifdef EYTA_DEBUG
void printKmersToFile (const Graph& graph, string filePath, string edgesPath);
//create the solid_kmers_binary_with_count file: a file with 2 columns, where the first is a k-mer encoded in 2-bit format and the second is the kmer abundance
void dump_solid_kmers(const gatb::core::debruijn::impl::Graph& graph, const string& kmer_counts_name);
#endif

#endif //KSGATB_DEBUG_H
