#include "map_reads.hpp"
#define NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME 1000
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>

struct MapReadsIteratorListener : public IteratorListener {
public:
    uint64_t &nbOfReadsProcessed;
    ISynchronizer* synchro;
    MapReadsIteratorListener(uint64_t &nbOfReadsProcessed, ISynchronizer* synchro) :
        nbOfReadsProcessed(nbOfReadsProcessed), synchro(synchro){}

    virtual void inc (u_int64_t ntasks_done) {
        // We lock the synchronizer
        synchro->lock ();

        nbOfReadsProcessed+=NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME;
        cerr << '\r' << nbOfReadsProcessed << " sequences mapped.";
        cerr.flush();

        // We unlock the synchronizer
        synchro->unlock ();
    }

    ~MapReadsIteratorListener(){}
};

//map the reads in the given file
void mapReadsMultithreaded(const vector <string> &allReadFilesNames, const string &descriptor, bool output,
                           const string &prefix, UnitigCounter &unitigCounter, int nbCores,
                           bool printReadInTheOutput, const SolidKmerOracle &solidKmerOracle,
                           int minSizeBlueUnitigs, int minSizeRedUnitigs, UnitigLinkingGraph* ulg, UnitigCounter* unitigCounterTranscripts, map<long, set<long>>* unitigId2ReadIndex) {

    //for each read file
    int fileIndex=0;
    for (const string &readFileName : allReadFilesNames) {
        // We declare an input Bank and use it locally
        IBank *inputBank = Bank::open(readFileName);
        LOCAL(inputBank);

        // Create the MappingOutputter
        MappingOutputter* mappingOutputter = NULL;
        if (output)
            mappingOutputter = new MappingOutputter(prefix+to_string(fileIndex));


        // We create a dispatcher configured for 'nbCores' cores and where each thread will map NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME reads at each time
        Dispatcher dispatcher(nbCores, NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME);

        //We create a progress iterator
        uint64_t nbOfReadsProcessed = 0;
        ISynchronizer* nbOfReadsProcessedSynchro = System::thread().newSynchronizer();

        //TODO: this new is deleted somewhere, it seems...
        MapReadsIteratorListener* mapReadsIteratorListener = new MapReadsIteratorListener(nbOfReadsProcessed, nbOfReadsProcessedSynchro);
        SubjectIterator<Sequence> it(inputBank->iterator(), NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME, mapReadsIteratorListener);

        cerr << "[Starting mapping " << descriptor << " - file: " << readFileName << " with " << nbCores << " cores... ]" << endl;
        cerr.flush();
        dispatcher.iterate(it,
                           MapAndPhaseSERead(fileIndex, mappingOutputter, prefix, nbOfReadsProcessed, nbOfReadsProcessedSynchro, unitigCounter,
                                       printReadInTheOutput, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs, ulg, unitigCounterTranscripts, unitigId2ReadIndex));
        cerr << endl << "[Mapping " << descriptor << " - file: " << readFileName << " finished!]" << endl;
        cerr.flush();

        if (output)
            delete mappingOutputter;

        fileIndex++;
    }
}



//map the reads in the given file in a PE-way
void mapPairedReadsMultithreaded(const vector <string> &allReadFilesNames, const string &descriptor, bool output,
                           const string &prefix, UnitigCounter &unitigCounter, int nbCores,
                           bool printReadInTheOutput, const SolidKmerOracle &solidKmerOracle,
                           int minSizeBlueUnitigs, int minSizeRedUnitigs, UnitigLinkingGraph* ulg, UnitigCounter* unitigCounterTranscripts) {

    //for each read file
    for (size_t fileIndex=0; fileIndex<allReadFilesNames.size(); fileIndex+=2) {
        string leftReadFileName=allReadFilesNames[fileIndex];
        string rightReadFileName=allReadFilesNames[fileIndex+1];

        // We declare the input Bank and use it locally
        IBank *leftInputBank = Bank::open(leftReadFileName);
        LOCAL(leftInputBank);
        IBank *rightInputBank = Bank::open(rightReadFileName);
        LOCAL(rightInputBank);

        //TODO: is this deleted somehow? scoped_ptr / unique_ptr bugs here... so I removed it..
        PairedIterator< Sequence, Sequence > * peIterator =
            new PairedIterator< Sequence, Sequence >(leftInputBank->iterator(), rightInputBank->iterator());

        // Create the MappingOutputter
        MappingOutputter* mappingOutputter = NULL;
        if (output)
            mappingOutputter = new MappingOutputter(prefix+to_string(fileIndex));


        // We create a dispatcher configured for 'nbCores' cores and where each thread will map NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME reads at each time
        Dispatcher dispatcher(nbCores, NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME);

        //We create a progress iterator
        uint64_t nbOfReadsProcessed = 0;
        ISynchronizer* nbOfReadsProcessedSynchro = System::thread().newSynchronizer();

        //TODO: this new is deleted somewhere, it seems...
        MapReadsIteratorListener* mapReadsIteratorListener = new MapReadsIteratorListener(nbOfReadsProcessed, nbOfReadsProcessedSynchro);
        SubjectIterator< pair<Sequence, Sequence> > it(peIterator, NB_OF_READS_EACH_THREAD_MAPS_AT_A_TIME, mapReadsIteratorListener);



        cerr << "[Starting mapping " << descriptor << " - left file: " << leftReadFileName << " , right file: "
             << rightReadFileName << " with " << nbCores << " cores... ]" << endl;
        cerr.flush();
        //dispatcher.iterate(it,
        dispatcher.iterate(peIterator,
                           MapAndPhasePERead(fileIndex, mappingOutputter, prefix, nbOfReadsProcessed, nbOfReadsProcessedSynchro, unitigCounter,
                                       printReadInTheOutput, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs, ulg, unitigCounterTranscripts, NULL));
        cerr << endl << "[Mapping " << descriptor << " - left file: " << leftReadFileName << " , right file: "
                     << rightReadFileName << " finished!]" << endl;
        cerr.flush();

        if (output)
            delete mappingOutputter;
    }
}

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
map_reads::map_reads ()  : Tool ("map_reads") //give a name to our tool
{
    populateParser(this);
}



void map_reads::checkParameters(){
    //check if we skip or not
    skip2 = getInput()->get(STR_PAR_SKIP2) != 0;

    //check if all files in STR_PAR_SE_SHORT_READS, STR_PAR_PE_SHORT_READS and STR_PAR_TRANSCRIPTOME exists
    build_dbg::checkReadFiles(this);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void map_reads::execute ()
{
    //map the transcripts to the graph
    checkParameters();
    if (skip2) return;

    string transcriptomeFile = getInput()->getStr(STR_PAR_TRANSCRIPTOME);
    string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
    string prefixMerged = prefix + "_merged";
    int nbCores = getInput()->getInt(STR_PAR_NB_CORES);
    int kmerSize = getInput()->getInt(STR_PAR_K);
    int minSizeBlueUnitigs = kmerSize + getInput()->getInt(STR_PAR_MIN_SIZE_BLUE_UNITIGS);
    int minSizeRedUnitigs = kmerSize + getInput()->getInt(STR_PAR_MIN_SIZE_RED_UNITIGS);

    SolidKmerOracle solidKmerOracle;
    {
        string solidKmerOracleFilename =  prefixMerged + ".sko";
        //load solidKmerOracle from the disk
        cerr << "Loading " << solidKmerOracleFilename << " ..." << endl;
        std::ifstream solidKmerOracleIF(solidKmerOracleFilename);
        boost::archive::text_iarchive textIArchive(solidKmerOracleIF);
        textIArchive >> solidKmerOracle;
        cerr << "Loading " << solidKmerOracleFilename << " ... - Done!" << endl;
    }


    //Do the Mapping and the Phasing, and colors the graph
    //Maps all the transcripts back to the graph
    //this will count how many transcripts map to an unitig
    UnitigCounter unitigCounterTranscripts(solidKmerOracle.getNbOfUnitigs(), 'B');

    //this will associate each blue unitig to the transcripts it is contained in
    map<long, set<long>> blueUnitigId2TranscriptIndex;

    //get all the transcriptomes
    vector <string> transcriptomeFileAsVector;
    transcriptomeFileAsVector.push_back(transcriptomeFile);

    //map all transcriptomes
    mapReadsMultithreaded(transcriptomeFileAsVector, "Transcriptome", true, prefix + "_merged.tr.mapping.",
                          unitigCounterTranscripts, nbCores, true, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs, NULL, NULL, &blueUnitigId2TranscriptIndex);

    #ifdef EYTA_DEBUG
    //cerr << "Test the mapping with scripts/check_if_mapping_is_ok.py"
    //exit(1);
    #endif


    //map the short reads file in order to build the ULG
    //this will count how many SRs map to an unitig
    UnitigCounter unitigCounterShortReads(solidKmerOracle.getNbOfUnitigs(), 'R');
    UnitigLinkingGraph ulg;


    //map the SINGLE END short reads unitigs
    //get all the short read files' name
    if (weHaveSEReads) {
        string SEReadsFile = getInput()->getStr(STR_PAR_SE_SHORT_READS);
        vector <string> SEShortReadFilesNames = getAllReadFilesNames(SEReadsFile);
        if (SEShortReadFilesNames.size()>0)
        {
            mapReadsMultithreaded(SEShortReadFilesNames, "SINGLE END short reads", false, prefix + "_merged.sesr.mapping.",
                                  unitigCounterShortReads, nbCores, false, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs,
                                  &ulg, &unitigCounterTranscripts, NULL);
        }
    }



    //map the PAIRED END short reads unitigs
    //get all the short read files' name
    if (weHavePEReads) {
        string PEReadsFile = getInput()->getStr(STR_PAR_PE_SHORT_READS);
        vector <string> PEShortReadFilesNames = getAllReadFilesNames(PEReadsFile);
        if (PEShortReadFilesNames.size()>0)
        {
            //TODO: put this back after debugging the PE reads
            /*
            mapPairedReadsMultithreaded(PEShortReadFilesNames, "PAIRED END short reads", false, "dummy",
                                        unitigCounterShortReads, nbCores, false, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs,
                                        &ulg, &unitigCounterTranscripts);
            */
            mapPairedReadsMultithreaded(PEShortReadFilesNames, "PAIRED END short reads", true, prefix + "_merged.pesr.mapping.",
                                        unitigCounterShortReads, nbCores, true, solidKmerOracle, minSizeBlueUnitigs, minSizeRedUnitigs,
                                        &ulg, &unitigCounterTranscripts);
        }
    }


    //serialize the ULG to a file
    {
        string ulgFilename=prefix+"_merged.ulg";
        cerr << "Saving " << ulgFilename << " ..." << endl;
        std::ofstream ulgOF(ulgFilename);
        boost::archive::text_oarchive textOArchive(ulgOF);
        textOArchive << ulg;
        cerr << "Saving " << ulgFilename << " ... - Done!" << endl;
    }

    //print the ulg to text files so we can see if it worked well or not
    //TODO: there is no need for this, this is just for DEBUG!!!
    //TODO: REMOVE ME!!!
    /*
    {
        ofstream ulgFile;
        openFileForWriting(prefix+"ulgEdges.txt", ulgFile);
        ulgFile << ulg.toStringEdges();
        ulgFile.close();
    }
    {
        ofstream ulgFile;
        openFileForWriting(prefix+"ulgHints.txt", ulgFile);
        ulgFile << ulg.toStringHints();
        ulgFile.close();
    }
     */
    //TODO: there is no need for this, this is just for DEBUG!!!
    //TODO: REMOVE ME!!!



    //output the node counts in a new file
    cerr << "Coloring and quantifying the graph..." << endl;
    {
        ifstream inputGraphFile;
        openFileForReading(prefix+"_merged.nodes", inputGraphFile);

        ofstream outputColoredGraphFile;
        openFileForWriting(prefix+"_merged.colored.nodes", outputColoredGraphFile);

        string line;
        while (getline(inputGraphFile, line)) {
            long nodeId; string nodeSeq;
            {
                stringstream ss;
                ss << line;
                ss >> nodeId >> nodeSeq;
            }
            outputColoredGraphFile << nodeId << "\t" <<
                solidKmerOracle.isTrustableUnitig(nodeId, &unitigCounterTranscripts, minSizeBlueUnitigs, minSizeRedUnitigs) << "\t" <<
                unitigCounterTranscripts.getNodeCount(nodeId) << "\t" <<
                unitigCounterShortReads.getNodeCount(nodeId) << "\t" <<
                nodeSeq << endl;
        }

        inputGraphFile.close();
        outputColoredGraphFile.close();


    }
    cerr << "Coloring and quantifying the graph... - Done!" << endl;

    //serialize _blue_trustable_nodes_to_transcripts
    cerr << "Saving " << (prefix+"_merged_blue_trustable_nodes_to_transcripts")  << "..." << endl;
    {
        ofstream blueTrustableNodesToTranscriptsFile(prefix + "_merged_blue_trustable_nodes_to_transcripts");
        boost::archive::text_oarchive textOArchive(blueTrustableNodesToTranscriptsFile);
        textOArchive << blueUnitigId2TranscriptIndex;
    }
    cerr << "Saving " << (prefix+"_merged_blue_trustable_nodes_to_transcripts")  << "... - Done!" << endl;
}
