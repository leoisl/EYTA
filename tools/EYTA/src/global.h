//
// Created by Leandro Ishi Soares de Lima on 19/07/16.
//

#ifndef KSGATB_GLOBAL_H
#define KSGATB_GLOBAL_H

#include <gatb/gatb_core.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "Utils.h"
#include "Strategy.h"
#include "SolidKmerOracle.h"
#include "Stats.h"
#define EYTA_DEBUG



//global vars
extern bool skip1;
extern bool skipTrGraph;
extern bool skipSrGraph;
extern bool skipMergedGraph;
extern bool skip2;
extern bool weHaveSEReads;
extern bool weHavePEReads;
extern Stats stats;


//parameters
extern const char* STR_PAR_SE_SHORT_READS;
extern const char* STR_PAR_PE_SHORT_READS;
extern const char* STR_PAR_K;
extern const char* STR_PAR_SHORT_READS_MIN_ABUNDANCE;
extern const char* STR_PAR_GRAPH_PREFIX;
extern const char* STR_PAR_NB_CORES;
extern const char* STR_PAR_SHORT_READS_RELATIVE_CUTOFF;
extern const char* STR_PAR_MIN_SIZE_BLUE_UNITIGS;
extern const char* STR_PAR_MIN_SIZE_RED_UNITIGS;
extern const char* STR_PAR_MAX_DEGREE_ULG;
extern const char* STR_PAR_TRANSCRIPTOME;
extern const char* STR_PAR_MAX_LENGTH_INTRON;
extern const char* STR_PAR_MAX_LENGTH_AS;
extern const char* STR_PAR_MIN_EDIT_DISTANCE;
extern const char* STR_PAR_SKIP1;
extern const char* STR_PAR_SKIP_SR_GRAPH;
extern const char* STR_PAR_SKIP_TR_GRAPH;
extern const char* STR_PAR_SKIP_MERGED_GRAPH;
extern const char* STR_PAR_SKIP2;
extern const char* STR_OUTPUT_CONTEXT;
extern const char* STR_AS_SPLICING_COMPLEXITY;
extern const char* STR_IR_SPLICING_COMPLEXITY;



void populateParser (Tool *tool);
#endif //KSGATB_GLOBAL_H
