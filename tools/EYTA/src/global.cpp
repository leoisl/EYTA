//
// Created by Leandro Ishi Soares de Lima on 04/02/17.
//

#include "global.h"

const char* STR_PAR_SE_SHORT_READS = "-sesr";
const char* STR_PAR_PE_SHORT_READS = "-pesr";
const char* STR_PAR_K = "-k";
const char* STR_PAR_SHORT_READS_MIN_ABUNDANCE = "-min_abundance";
const char* STR_PAR_GRAPH_PREFIX = "-prefix";
const char* STR_PAR_NB_CORES = "-nb-cores";
const char* STR_PAR_SHORT_READS_RELATIVE_CUTOFF = "-rel_cutoff";
const char* STR_PAR_MIN_SIZE_BLUE_UNITIGS = "-minSizeBlueUnitigs";
const char* STR_PAR_MIN_SIZE_RED_UNITIGS = "-minSizeRedUnitigs";
const char* STR_PAR_TRANSCRIPTOME = "-t";
const char* STR_PAR_MAX_LENGTH_INTRON = "-max_length_intron";
const char* STR_PAR_MAX_LENGTH_AS = "-max_length_as";
const char* STR_PAR_MIN_EDIT_DISTANCE = "-min_edit_distance";
const char* STR_PAR_SKIP1 = "-skip1";
const char* STR_PAR_SKIP_SR_GRAPH = "-skipSR";
const char* STR_PAR_SKIP_TR_GRAPH = "-skipTR";
const char* STR_PAR_SKIP_MERGED_GRAPH = "-skipMerged";
const char* STR_PAR_SKIP2 = "-skip2";
const char* STR_OUTPUT_CONTEXT = "-output_context";
const char* STR_AS_SPLICING_COMPLEXITY = "-splicing_complexity_AS";
const char* STR_IR_SPLICING_COMPLEXITY = "-splicing_complexity_IR";


bool skip1 = false;
bool skipTrGraph = false;
bool skipSrGraph = false;
bool skipMergedGraph = false;
bool skip2 = false;
bool weHaveSEReads = false;
bool weHavePEReads = false;
Stats stats;


void populateParser (Tool *tool) {
  //TODO: fix -nb-cores=1 by default

  // We add some custom arguments for command line interface
  tool->getParser()->push_front (new OptionNoParam (STR_PAR_SKIP1, "Skips Step 1 entirely, running only Step 2. Assumes that Step 1 was correctly run.",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_PAR_SKIP_SR_GRAPH, "Skips building the short reads graph in Step 1. Uses pre-built short reads graph.",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_PAR_SKIP_TR_GRAPH, "Skips building the transcriptome graph in Step 1. Uses pre-built transcriptome graph.",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_PAR_SKIP_MERGED_GRAPH, "Skips building the merged graph in Step 1. Uses pre-built merged graph.",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_PAR_SKIP2, "Skips Step 2 entirely (transcriptome mapping and ULG construction). Uses pre-built transcriptome mapping and ULG.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_GRAPH_PREFIX, "Prefix of the name of the built files related to the graph",  false, "graph"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_SHORT_READS_MIN_ABUNDANCE, "K-mers present strictly less than this"
      "number of times in the dataset will be discarded. This is applied only on the short reads graph.",  false, "2"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_SHORT_READS_RELATIVE_CUTOFF, "A threshold. Edges which counts are relatively smaller than this threshold are removed. This is applied only on the short reads graph only.",  false, "0.02"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_MIN_SIZE_BLUE_UNITIGS, "An integer, min_blue. The minimum size of a blue unitig to be used in assembly will be k+min_blue",  false, "5"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_MIN_SIZE_RED_UNITIGS, "An integer, min_red. The minimum size of a red-only unitig to be used in assembly will be k+min_red",  false, "15"));
  tool->getParser()->push_front (new OptionNoParam (STR_OUTPUT_CONTEXT, "If we should output the context or not of the flanking unitigs", false));

  //I've set this maximum length to 11100 (this comes from the fact that < 10% of human introns are more than 11000 bp in length
    //see https://www.researchgate.net/publication/8491627_Distributions_of_exons_and_introns_in_the_human_genome )
  //And due to the small illumina error rate and some insertion variants, 11100 seemed a nice guess to try to capture introns also
  //TODO: check this!!
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_MAX_LENGTH_INTRON, "Maximum length of the searched introns",  false, "11100"));
  //TODO: check this!!

  //This concerns the maximum length of the searched AS events excluding introns
  //Based on https://www.researchgate.net/publication/8491627_Distributions_of_exons_and_introns_in_the_human_genome, we can see that the majority of human exons have <= 600 bp in length
  //To be able to capture multiple exon skipping events, up to 3, we have set this to 2000
  //TODO: check this!!
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_MAX_LENGTH_AS, "Maximum length of the searched alternative splicing events (except introns)",  false, "2000"));
  //TODO: check this!!

  tool->getParser()->push_front (new OptionOneParam (STR_IR_SPLICING_COMPLEXITY, "The maximum number of alternative paths to find between two constitutive exons - when searching for introns - WARNING as introns usually traverses a lot of repeats, increasing this value can be dangerous - 1 means that we search for an unique path that makes the most sense.",  false, "1"));
  tool->getParser()->push_front (new OptionOneParam (STR_AS_SPLICING_COMPLEXITY, "The maximum number of alternative paths to find between two constitutive exons - when searching for AS except intron (exon skipping, alternative acceptor, alternative donor, etc...)",  false, "10"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_MIN_EDIT_DISTANCE, "Minimum relative edit distance threshold to consider that two sequences represent different splicing events.", false, "0.2"));
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_K, "K-mer size",  false, "31"));

  tool->getParser()->push_front (new OptionOneParam (STR_PAR_SE_SHORT_READS, "A text file containing one SINGLE END short read file at each line. Can be used together with parameter -pesr if you have single and paired end reads. At least -sesr or -pesr must be specified.", false, ""));

  //TODO: add it back after!!!
  //removing PE-reads from the input
  //tool->getParser()->push_front (new OptionOneParam (STR_PAR_PE_SHORT_READS, "A text file containing PAIRED END short read files. The 1st file (left) is paired with the 2nd (right), the 3rd file (left) is paired with the 4th (right), and so on... Can be used together with parameter -sesr if you have single and paired end reads. At least -sesr or -pesr must be specified.",  false, ""));
  //removing PE-reads from the input
  //TODO: add it back after!!!
  
  tool->getParser()->push_front (new OptionOneParam (STR_PAR_TRANSCRIPTOME, "A path to a FASTA file representing a transcriptome.",  true));
}

