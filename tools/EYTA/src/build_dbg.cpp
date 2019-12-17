//! [snippet1]

#include "build_dbg.hpp"
#include "global.h"
#include "GraphOutput.h"
#include "Tests.h"
#include <boost/archive/text_oarchive.hpp>

#ifdef EYTA_DEBUG
#include "debug.h"
#endif

using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
build_dbg::build_dbg ()  : Tool ("build_dbg") //give a name to our tool
{
    populateParser(this);
}

void buildSequence (
        const Graph& graph,
        const Node& startingNode,
        size_t length,
        size_t nbContigs,
        const string& consensusRight,
        const string& consensusLeft,
        Sequence& seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();

    /** We set the sequence comment. */
    stringstream ss1;
    ss1 << nbContigs << "__len__" << length << " ";
    seq._comment = ss1.str();

    /** We set the data length. */
    seq.getData().resize (length);

    //We fill the data
    string finalSequence = consensusLeft + graph.toString (startingNode) + consensusRight;
    for (size_t i=0;i<finalSequence.size();i++)
        data[i] = finalSequence[i];
}

void construct_linear_seqs (const gatb::core::debruijn::impl::Graph& graph, const string& prefix, SolidKmerOracle* solidKmerOracle = NULL)
{
    using namespace gatb::core::debruijn::impl;
    using namespace gatb::core::tools::misc::impl;

    string linear_seqs_name = prefix+".unitigs";
    string unitigCoverageFilename = linear_seqs_name+".coverage";
    string unitigSizeFilename = linear_seqs_name+".size";

    IBank* outputBank = new BankFasta (linear_seqs_name.c_str());
    LOCAL (outputBank);
    ofstream unitigCoverageFile(unitigCoverageFilename);

    // We create a Terminator object - this will mark the nodes that are already built
    MPHFTerminator terminator (graph);

    // We create a BranchingTerminator object - this will mark the nodes where to stop the traversal
    BranchingTerminator branchingTerminator(graph);

    // We create a Traversal instance to traverse unitigs
    Traversal* traversal = Traversal::create (TRAVERSAL_UNITIG, graph, branchingTerminator);
    LOCAL (traversal);

    Path consensusRight;
    Path consensusLeft;
    Sequence seq (Data::ASCII);
    u_int64_t nbContigs = 0;
    BankFasta::setDataLineSize(0);

    //We loop through the nodes and build the unitigs
    ProgressGraphIterator<Node, ProgressTimerAndSystem> it (graph.iterator(), "Graph: building unitigs");
    for (it.first(); !it.isDone(); it.next()) {
        auto &startingNode = it.item();

        if (terminator.is_marked(startingNode))
            continue;

        auto reversedNode = graph.reverse(startingNode);
        int lenRight = traversal->traverse (startingNode, DIR_OUTCOMING, consensusRight);
        int lenLeft = traversal->traverse (reversedNode, DIR_OUTCOMING, consensusLeft);
        int lenTotal = graph.getKmerSize() + lenRight + lenLeft;

        //mark the traversed nodes
        terminator.mark(startingNode);
        auto currentNode = startingNode;
        for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });
        currentNode = reversedNode;
        for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });
        auto leftestNode = graph.reverse(currentNode);

        // We get the unitig strings
        string consensusLeftStr;
        {
            stringstream ss;
            ss << consensusLeft;
            consensusLeftStr=reverse_complement(ss.str());
        }

        string consensusRightStr;
        {
            stringstream ss;
            ss << consensusRight;
            consensusRightStr=ss.str();
        }


        /** We create the contig sequence. */
        buildSequence(graph, startingNode, lenTotal, nbContigs, consensusRightStr, consensusLeftStr, seq);

        if (solidKmerOracle) { //if it was given, use it
            //we add it to the oracle
            solidKmerOracle->addUnitigSequence(seq.toString(), leftestNode);
        }

        /** We add the sequence into the output bank. */
        outputBank->insert (seq);


        //We calculate the coverage of this unitig, i.e. the average k-mer count
        vector<CountNumber> kmerCoverages;
        kmerCoverages.push_back(graph.queryAbundance(startingNode));
        currentNode = startingNode;

        for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            kmerCoverages.push_back(graph.queryAbundance(currentNode));
        });
        currentNode = graph.reverse(startingNode);
        for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            kmerCoverages.push_back(graph.queryAbundance(currentNode));
        });

        //compute the average
        double coverage=0;
        for(const auto &kmerCoverage : kmerCoverages)
            coverage+=kmerCoverage;
        coverage=coverage/kmerCoverages.size();

        //output the coverage
        //TODO: binary writing
        unitigCoverageFile << coverage << endl;

        //increase the number of contigs
        nbContigs += 1;
    }

    //print the nb of contigs to unitigSizeFilename file
    ofstream unitigSizeFile;
    openFileForWriting(unitigSizeFilename, unitigSizeFile);
    unitigSizeFile << nbContigs;

    //close all files
    outputBank->flush ();
    unitigSizeFile.close();
    unitigCoverageFile.close();
}



class EdgeConstructionVisitor : public boost::static_visitor<>    {
private:
    const string& linear_seqs_name;

public:
    EdgeConstructionVisitor (const string &linear_seqs_name) : linear_seqs_name(linear_seqs_name) {}
    template<size_t span>
    void operator() (GraphOutput<span>& graphOutput) const
    {
        graphOutput.open();
        graphOutput.load_nodes_extremities(linear_seqs_name);
        graphOutput.construct_graph(linear_seqs_name);
        graphOutput.close();
    }
};

//check if all files in STR_PAR_SE_SHORT_READS, STR_PAR_PE_SHORT_READS and STR_PAR_TRANSCRIPTOME exists
void build_dbg::checkReadFiles(Tool* tool) {
    //check the read files
    if (tool->getInput()->get(STR_PAR_SE_SHORT_READS)!=0)
        weHaveSEReads=true;
    if (tool->getInput()->get(STR_PAR_PE_SHORT_READS)!=0)
        weHavePEReads=true;

    if (!weHaveSEReads && !weHavePEReads) {
        stringstream ss;
        ss << "[PARAMETER ERROR] You should provide at least " << STR_PAR_SE_SHORT_READS << " or " << STR_PAR_PE_SHORT_READS << endl;
        fatalError(ss.str());
    }

    //check if all files in STR_PAR_SE_SHORT_READS, STR_PAR_PE_SHORT_READS and STR_PAR_TRANSCRIPTOME exists
    if (weHaveSEReads) {
        string SEReads = tool->getInput()->getStr(STR_PAR_SE_SHORT_READS);
        checkIfReadFileIsFine(SEReads);
    }
    if (weHavePEReads) {
        string PEReads = tool->getInput()->getStr(STR_PAR_PE_SHORT_READS);
        checkIfReadFileIsFine(PEReads);
    }
    {
        string transcriptsFile = tool->getInput()->getStr(STR_PAR_TRANSCRIPTOME);
        checkIfFileExists(transcriptsFile);
    }
}

void build_dbg::checkParameters() {
    //check if we skip or not
    skip1 = getInput()->get(STR_PAR_SKIP1) != 0;
    skipSrGraph = getInput()->get(STR_PAR_SKIP_SR_GRAPH) != 0;
    skipTrGraph = getInput()->get(STR_PAR_SKIP_TR_GRAPH) != 0;
    skipMergedGraph = getInput()->get(STR_PAR_SKIP_MERGED_GRAPH) != 0;

    //check if all files in STR_PAR_SE_SHORT_READS, STR_PAR_PE_SHORT_READS and STR_PAR_TRANSCRIPTOME exists
    checkReadFiles(this);
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void build_dbg::execute ()
{
    //check and get parameters
    checkParameters();

    //skip step 1 completely
    if (skip1) return;

    //read parameters
    string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
    string prefixSR = prefix + "_sr";
    string prefixTr = prefix + "_tr";
    string prefixMerged = prefix + "_merged";

    //build the file with all the reads
    string allReadsFile;
    if (weHaveSEReads && weHavePEReads) {
        //we have both files, create a file containing all
        string SEReadsFile = getInput()->getStr(STR_PAR_SE_SHORT_READS);
        string PEReadsFile = getInput()->getStr(STR_PAR_PE_SHORT_READS);
        allReadsFile = prefix+"_allReadsFile";
        concatenate2FilesIntoA3rd(SEReadsFile, PEReadsFile, allReadsFile);
    }else if (weHaveSEReads) {
        //we have just SE reads file
        allReadsFile = getInput()->getStr(STR_PAR_SE_SHORT_READS);
    }else {
        //we have just PE reads file
        allReadsFile = getInput()->getStr(STR_PAR_PE_SHORT_READS);
    }

    int kmerSize = getInput()->getInt(STR_PAR_K);
    int minAbundance = getInput()->getInt(STR_PAR_SHORT_READS_MIN_ABUNDANCE);
    int nbCores = getInput()->getInt(STR_PAR_NB_CORES);
    double relativeCutoff = 0;
    if (getInput()->get(STR_PAR_SHORT_READS_RELATIVE_CUTOFF))
        relativeCutoff = getInput()->getDouble(STR_PAR_SHORT_READS_RELATIVE_CUTOFF);


    /*
    [graph options]

   [kmer count options]
          -in                             (1 arg) :    reads file
          -kmer-size                      (1 arg) :    size of a kmer  [default '31']
          -abundance-min                  (1 arg) :    min abundance threshold for solid kmers  [default '2']
          -abundance-max                  (1 arg) :    max abundance threshold for solid kmers  [default '2147483647']
          -abundance-min-threshold        (1 arg) :    min abundance hard threshold (only used when min abundance is "auto")  [default '2']
          -histo-max                      (1 arg) :    max number of values in kmers histogram  [default '10000']
          -solidity-kind                  (1 arg) :    way to compute counts of several files (sum, min, max, one, all, custom)  [default 'sum']
          -solidity-custom                (1 arg) :    when solidity-kind is cutom, specifies list of files where kmer must be present  [default '']
          -max-memory                     (1 arg) :    max memory (in MBytes)  [default '5000']
          -max-disk                       (1 arg) :    max disk   (in MBytes)  [default '0']
          -solid-kmers-out                (1 arg) :    output file for solid kmers (only when constructing a graph)  [default '']
          -out                            (1 arg) :    output file  [default '']
          -out-dir                        (1 arg) :    output directory  [default '.']
          -out-tmp                        (1 arg) :    output directory for temporary files  [default '.']
          -out-compress                   (1 arg) :    h5 compression level (0:none, 9:best)  [default '0']
          -storage-type                   (1 arg) :    storage type of kmer counts ('hdf5' or 'file')  [default 'hdf5']

      [kmer count, algorithmic options options]
             -minimizer-type   (1 arg) :    minimizer type (0=lexi, 1=freq)  [default '0']
             -minimizer-size   (1 arg) :    size of a minimizer  [default '10']
             -repartition-type (1 arg) :    minimizer repartition (0=unordered, 1=ordered)  [default '0']

   [bloom options]
          -bloom        (1 arg) :    bloom type ('basic', 'cache', 'neighbor')  [default 'neighbor']
          -debloom      (1 arg) :    debloom type ('none', 'original' or 'cascading')  [default 'cascading']
          -debloom-impl (1 arg) :    debloom impl ('basic', 'minimizer')  [default 'minimizer']

   [branching options]
          -branching-nodes (1 arg) :    branching type ('none' or 'stored')  [default 'stored']
          -topology-stats  (1 arg) :    topological information level (0 for none)  [default '0']

   [general options]
          -config-only       (0 arg) :    dump config only
          -nb-cores          (1 arg) :    number of cores  [default '0']
          -verbose           (1 arg) :    verbosity level  [default '1']
          -integer-precision (1 arg) :    integers precision (0 for optimized value)  [default '0']
     */

    //Builds the SR DBG using GATB
    if(!skipSrGraph)
    {
        //TODO: by using create() and assigning to a Graph object, the copy constructor does a shallow or deep copy??
        //TODO: bug - if you run with -abundance-min 100, where we have 0 solid kmers, this bugs
        cerr << "Building Short Reads graph..." << endl;
        Graph srGraph = gatb::core::debruijn::impl::Graph::create(
            "-in %s -kmer-size %d -abundance-min %d -out %s -nb-cores %d",
            allReadsFile.c_str(), kmerSize, minAbundance, prefixSR.c_str(), nbCores);

        //computes the edge list for each node
        srGraph.precomputeAdjacency(nbCores, true);


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //DEBUG
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        #ifdef EYTA_DEBUG
        //test if the error removal is ok
        //EYTATester::relativeErrorRemovalTest(srGraph, nbCores, prefixSR, relativeCutoff);
        #endif
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //DEBUG
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //do the error removal by removing edges having relatively small count
        if (relativeCutoff > 0)
            srGraph.relativeErrorRemoval(relativeCutoff, prefixSR, nbCores);
        else {
            //create empty prefixSR + ".removed_edges" file
            ofstream dummyOfstream;
            openFileForWriting(prefixSR + ".removed_edges", dummyOfstream);
            dummyOfstream.close();
        }

        // building the unitigs
        construct_linear_seqs(srGraph, prefixSR);

        //cleanup
        boost::filesystem::remove(prefixSR + ".h5");

        cerr << "Building Short Reads graph - Done!" << endl;
    }else {
        cerr << "[SKIPPED] Building Short Reads graph" << endl;
    }

    ///Builds the transcriptome DBG using GATB
    if(!skipTrGraph)
    {
        //build the transcriptome graph
        cerr << "Building Transcriptome graph..." << endl;
        string transcriptsFile = getInput()->getStr(STR_PAR_TRANSCRIPTOME);
        Graph trGraph = gatb::core::debruijn::impl::Graph::create(
            "-in %s -kmer-size %d -abundance-min 0 -out %s -nb-cores %d",
            transcriptsFile.c_str(), kmerSize, prefixTr.c_str(), nbCores);

        //remove erroneously removed edges
        cerr << "Removing from " << prefixSR << ".removed_edges the edges that are correct according to the transcriptome..." << endl;

        //read the edges removed when building the short read graph into a set
        set<string> SRGraphRemovedEdges = readRemovedEdgesIntoASetOfString(prefixSR + ".removed_edges");

        //get all the edges in the transcriptomic graph that are in branching nodes
        set<string> edgesInTheTranscriptomicGraph;
        {
            ProgressGraphIteratorTemplate<BranchingNode_t<Node>, ProgressTimerAndSystem> itNode(
                trGraph.iteratorBranching(),
                "Getting all edges from the transcriptomic graph in branching nodes");
            for (itNode.first(); !itNode.isDone(); itNode.next()) {
                Node &node = itNode.item();

                for (Direction dir = DIR_OUTCOMING;
                     dir < DIR_END; dir = (Direction) ((int) dir + 1)) // in both directions
                {
                    auto neighbors = trGraph.neighborsEdge(node, dir);
                    for (unsigned int i = 0; i < neighbors.size(); i++) {
                        Edge edge = trGraph.getTheEquivalentOutgoingEdge(neighbors[i]);
                        edgesInTheTranscriptomicGraph.insert(trGraph.toString(edge));
                    }
                }
            }
        }

        //removes from SRGraphRemovedEdges the edges in edgesInTheTranscriptomicGraph
        cerr << "Removing edges..." << endl;
        //this contains the correctly removed edges
        vector<string> correctlyRemovedEdgesVector(SRGraphRemovedEdges.size());
        auto itSetDifference = set_difference(SRGraphRemovedEdges.begin(), SRGraphRemovedEdges.end(),
                       edgesInTheTranscriptomicGraph.begin(), edgesInTheTranscriptomicGraph.end(),
                       correctlyRemovedEdgesVector.begin());
        correctlyRemovedEdgesVector.resize(itSetDifference-correctlyRemovedEdgesVector.begin());

        cerr << (SRGraphRemovedEdges.size() - correctlyRemovedEdgesVector.size()) << " edges added back to the graph due to being in the transcriptome!" << endl;

        //print the edges to a file
        ofstream correctlyRemovedEdgesFile;
        openFileForWriting(prefixSR+"_"+prefixTr+".removed_edges", correctlyRemovedEdgesFile);
        for (const auto &edge : correctlyRemovedEdgesVector)
            correctlyRemovedEdgesFile << edge << endl;
        correctlyRemovedEdgesFile.close();
        cerr << "Removing edges... - Done!" << endl;
        cerr << "Removing from " << prefixSR << ".removed_edges the edges that are correct according to the transcriptome... - Done!" << endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //DEBUG
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        #ifdef EYTA_DEBUG
        //EYTATester::checkIfAddingBackTranscriptomicEdgeIsCorrect(SRGraphRemovedEdges, edgesInTheTranscriptomicGraph, correctlyRemovedEdges);
        #endif
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //DEBUG
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Finding the unitigs
        construct_linear_seqs(trGraph, prefixTr);

        //cleanup
        boost::filesystem::remove(prefixTr + ".h5");

        cerr << "Building Transcriptome graph - Done!" << endl;
    }else {
        cerr << "[SKIPPED] Building Transcriptome graph" << endl;
    }


    //Builds the DBG from both the unitigs file created by the two previous step
    if (skipMergedGraph == false)
    {
        cerr << "Building Merged graph..." << endl;

        //Builds a file containing both unitigs file to give to GATB
        string bothUnitigsFilepath = prefixMerged+".merged_unitigs_file";
        {
            ofstream bothUnitigsFile(bothUnitigsFilepath);
            bothUnitigsFile << prefixSR << ".unitigs" << endl
                            << prefixTr << ".unitigs" << endl;
            bothUnitigsFile.close();
        }

        //TODO: has to be created as pointer since SolidKmerOracle will reference it
        //SolidKmerOracle is responsible for deleting it... take care of this
        Graph* graph = gatb::core::debruijn::impl::Graph::createAsPointer(
            "-in %s -kmer-size %d -abundance-min 0 -out %s -nb-cores %d",
            bothUnitigsFilepath.c_str(), kmerSize, prefixMerged.c_str(), nbCores);

        //computes the edge list for each node
        graph->precomputeAdjacency(nbCores, true);

        //removes from the graph the edges that should be removed
        #ifdef EYTA_DEBUG
        //EYTATester::checkRemoveEdgesFromFile(*graph, prefixSR, prefixTr, nbCores, correctlyRemovedEdges);
        #endif
        graph->removeEdgesFromFile(prefixSR+"_"+prefixTr+".removed_edges", nbCores);

        //compute the unitigs of the merged graph
        string linear_seqs_name = prefixMerged+".unitigs";
        SolidKmerOracle solidKmerOracle(graph, prefixMerged); //also computes the solidKmerOracle
        construct_linear_seqs (*graph, prefixMerged, &solidKmerOracle);

        //builds and outputs .nodes and .edges.dbg files
        //TODO: go until GraphOutput<KMER_SPAN(7)>?
        typedef boost::variant <
            GraphOutput<KMER_SPAN(0)>,
            GraphOutput<KMER_SPAN(1)>
            /* TODO: ADD THIS BACK WHEN CHANGING THE KMER SIZE IN THE CMAKE_LIST */
            /*,
            GraphOutput<KMER_SPAN(2)>,
            GraphOutput<KMER_SPAN(3)>*/
        >  GraphOutputVariant;

        GraphOutputVariant graphOutput;
        //TODO: check this
        if (kmerSize < KMER_SPAN(0))  {  graphOutput = GraphOutput<KMER_SPAN(0)>(graph, prefixMerged); }
        else if (kmerSize < KMER_SPAN(1))  {  graphOutput = GraphOutput<KMER_SPAN(1)>(graph, prefixMerged); }
        /* TODO: ADD THIS BACK WHEN CHANGING THE KMER SIZE IN THE CMAKE_LIST
        else if (kmerSize < KMER_SPAN(2))  {  graphOutput = GraphOutput<KMER_SPAN(2)>(graph, prefixMerged); }
        else if (kmerSize < KMER_SPAN(3))  {  graphOutput = GraphOutput<KMER_SPAN(3)>(graph, prefixMerged); }
         */
        else { throw gatb::core::system::Exception ("Graph failure because of unhandled kmer size %d", kmerSize); }
        boost::apply_visitor (EdgeConstructionVisitor(linear_seqs_name),  graphOutput);

        #ifdef EYTA_DEBUG
        //EYTATester::testSolidKmerOracle(solidKmerOracle);
        #endif

        //save solidKmerOracle to the disk
        cerr << "Saving " << prefixMerged << ".sko ..." << endl;
        std::ofstream solidKmerOracleOF(prefixMerged+".sko");
        boost::archive::text_oarchive textOArchive(solidKmerOracleOF);
        textOArchive << solidKmerOracle;
        cerr << "Saving " << prefixMerged << ".sko ... - Done!" << endl;

        //save disk space
        boost::filesystem::remove(linear_seqs_name);
        cerr << "Building Merged graph - Done!" << endl;
    }else {
        cerr << "[SKIPPED] Building Merged graph" << endl;
    }
}
