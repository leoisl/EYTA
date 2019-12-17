//
// Created by Leandro Ishi Soares de Lima on 24/02/17.
//

#include "Tests.h"

/*
void kisSpliceErrorRemovalTest() {
  vector<IBank*> ibanks;
  vector<set<string> > correctEdges;

  // We get a handle on a fake bank made of 3 sequences.
  ibanks.push_back(new BankStrings (
      "AAAAAAAAAAAAAAAAAAAC",
      "AATGTGTCTGTCGATTGTCCT",
      "CATGTGTCTGTCGATTGTCCG",
      "AATGTGTCTGTCGATTGTCC",
      "AATGTGTCTGTCGATTGTCC",
      "AATGTGTCTGTCGATTGTCC",
      NULL
  ));
  // We create the graph with the bank and other options
  Graph graph = Graph::create(ibanks[0], "-kmer-size 19  -abundance-min 1 -mphf emphf -out testgraph");

  //computes the edge list for each node
  graph.precomputeAdjacency(1, true);

  auto it = graph.iterator();
  for (it.first(); !it.isDone(); it.next()) {
    Node node = it.item();
    for (Direction dir = DIR_OUTCOMING; dir < DIR_END; dir = (Direction) ((int) dir + 1)) // in both directions
    {
      GraphTemplate<Node, Edge, GraphDataVariant>::Vector<Edge> neighbors = graph.neighborsEdge(node, dir);
      for (unsigned int i = 0; i < neighbors.size(); i++) {
        cout << "Edge: " << graph.toString(neighbors[i]) << endl;
        cout << "getTheEquivalentOutgoingEdge: " << graph.toString(graph.getTheEquivalentOutgoingEdge(neighbors[i])) << endl;
      }
    }
  }
  exit(1);




  // We get a handle on a fake bank made of 3 sequences.
  ibanks.push_back(new BankStrings (
      "AATGTGTCTGTCGATTGTCCT",
      "CATGTGTCTGTCGATTGTCCG",
      "ATGTGTCTGTCGATTGTCCT",
      "ATGTGTCTGTCGATTGTCCT",
      "ATGTGTCTGTCGATTGTCCT",
      NULL
  ));
  {
    set<string> correctEdgesSet;
    correctEdgesSet.insert("4F 0F");
    correctEdgesSet.insert("0F 1R");
    correctEdgesSet.insert("0R 4R");
    correctEdgesSet.insert("0R 3R");
    correctEdgesSet.insert("1F 0R");
    correctEdgesSet.insert("3F 0F");
    correctEdges.push_back(correctEdgesSet);
  }

  ibanks.push_back(new BankStrings (
      "AATGTGTCTGTCGATTGTCCT",
      "CATGTGTCTGTCGATTGTCCG",
      "AATGTGTCTGTCGATTGTCC",
      "AATGTGTCTGTCGATTGTCC",
      "AATGTGTCTGTCGATTGTCC",
      NULL
  ));
  {
    set<string> correctEdgesSet;
    correctEdgesSet.insert("4F 0F");
    correctEdgesSet.insert("0F 1R");
    correctEdgesSet.insert("0R 4R");
    correctEdgesSet.insert("0R 3R");
    correctEdgesSet.insert("1F 0R");
    correctEdgesSet.insert("3F 0F");
    correctEdges.push_back(correctEdgesSet);
  }

  int i=0;
  for(auto bank : ibanks) {
    try {
      // We create the graph with the bank and other options
      Graph graph = Graph::create(bank, "-kmer-size 19  -abundance-min 1 -mphf emphf -out testgraph");

      //computes the edge list for each node
      graph.precomputeAdjacency(1, true);


      cout << "Writing edges before error removal: " << endl;
      set<string> allEdges;
      auto it = graph.iterator();
      for (it.first(); !it.isDone(); it.next()) {
        cout << "Node: " << graph.toString(it.item()) << " - index: " << graph.nodeMPHFIndex(it.item()) << " abundance: " << graph.queryAbundance(it.item()) << endl;
        graph.neighborsEdge(it.item().kmer).iterate([&](const Edge &e) {
            //cout << "Edge: " << graph.toString(e) << endl;
            Node from = e.from;
            Node to = e.to;
            cout << "Edge: " << graph.nodeMPHFIndex(from) << (from.strand==STRAND_FORWARD?'F':'R') << " " << graph.nodeMPHFIndex(to) << (to.strand==STRAND_FORWARD?'F':'R') << endl;
        });
      }
      cout << endl;

      //do the error removal by removing edges from the precomputed adjacency list
      graph.relativeErrorRemoval(0.5, 1);

      cout << "Writing edges after error removal: " << endl;
      it = graph.iterator();
      set<string> edgesWithoutErrors;

      for (it.first(); !it.isDone(); it.next()) {
        cout << "Node: " << graph.toString(it.item()) << " - index: " << graph.nodeMPHFIndex(it.item()) << " abundance: " << graph.queryAbundance(it.item()) << endl;
        graph.neighborsEdge(it.item().kmer).iterate([&](const Edge &e) {
            //cout << "Edge: " << graph.toString(e) << endl;
            Node from = e.from;
            Node to = e.to;

            string edgeAsStr="";
            {
              stringstream ss;
              ss << graph.nodeMPHFIndex(from) << (from.strand==STRAND_FORWARD?'F':'R') << " " << graph.nodeMPHFIndex(to) << (to.strand==STRAND_FORWARD?'F':'R');
              edgeAsStr = ss.str();
            }

            cout << "Edge: "  << edgeAsStr << endl;
            edgesWithoutErrors.insert(edgeAsStr);
        });
      }

      if (edgesWithoutErrors!=correctEdges[i]) {
        cout << "BUG ON " << i << endl;
        exit(1);
      } else {
        cout << "TEST OK!!!" << endl;
      }

    }
    catch (Exception &e) {
      std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
    i++;
  }
}
 */


void relativeErrorRemovalTest(Graph &graph, int nbCores, const string &prefix, double relativeCutoff) {
  cerr << "Starting relativeErrorRemovalTest()..." << endl;

  //get all the edges that were in the graph
  set<string> edgesBefore;
  auto it = graph.iterator();
  for (it.first(); !it.isDone(); it.next()) {
    Node node = it.item();
    for (Direction dir = DIR_OUTCOMING; dir < DIR_END; dir = (Direction) ((int) dir + 1)) // in both directions
    {
      auto neighbors = graph.neighborsEdge(node, dir);
      for (unsigned int i = 0; i < neighbors.size(); i++) {
        edgesBefore.insert(graph.toString(graph.getTheEquivalentOutgoingEdge(neighbors[i])));
      }
    }
  }


  //do the error removal
  if (relativeCutoff > 0)
    graph.relativeErrorRemoval(relativeCutoff, prefix, nbCores);


  //get all the edges after the error removal
  set<string> edgesAfter;
  it = graph.iterator();
  for (it.first(); !it.isDone(); it.next()) {
    Node node = it.item();
    for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) ) // in both directions
    {
      auto neighbors = graph.neighborsEdge(node, dir);
      for (unsigned int i = 0; i < neighbors.size(); i++) {
        edgesAfter.insert(graph.toString(graph.getTheEquivalentOutgoingEdge(neighbors[i])));
      }
    }
  }


  //get the file with the removed edges
  ifstream removedEdgesFile(prefix+".removed_edges");
  set<string> removedEdges;
  string line;
  while(getline(removedEdgesFile, line))
    removedEdges.insert(line);
  removedEdgesFile.close();


  //check if they are equal
  for(string removedEdge : removedEdges) {
    cerr << "[relativeErrorRemovalTest] Testing " << removedEdge << endl;
    if (edgesBefore.find(removedEdge) == edgesBefore.end()) {
      cerr << "[ERROR]: " << removedEdge << " NOT FOUND in edges before ER." << endl;
      exit(1);
    }
    if (edgesAfter.find(removedEdge) != edgesAfter.end()) {
      cerr << "[ERROR]: " << removedEdge << " FOUND in edges after ER, but it should not (it was removed)." << endl;
      exit(1);
    }
  }

  cerr << "relativeErrorRemovalTest() done, no errors!" << endl;
  exit(1);
}



/*
void printAllOutgoingEdgeToFile (const string &filename, Graph* graph) {
  ofstream file(filename);
  auto it = graph->iterator();
  for (it.first(); !it.isDone(); it.next()) {
    Node node = it.item();
    for (Direction dir = DIR_OUTCOMING; dir < DIR_END; dir = (Direction) ((int) dir + 1)) // in both directions
    {
      auto neighbors = graph->neighborsEdge(node, dir);
      for (unsigned int i = 0; i < neighbors.size(); i++) {
        file << graph->toString(graph->getTheEquivalentOutgoingEdge(neighbors[i])) << endl;
      }
    }
  }
  file.close();
}
 */


void EYTATester::testSolidKmerOracle(SolidKmerOracle* solidKmerOracle) {
  cout << "EYTATester::testSolidKmerOracle() begin test..." << endl;
  cout.flush();

  GraphIterator<Node> it = solidKmerOracle->graph->iterator();
  int kmersWithoutUnitigs=0;
  for (it.first(); !it.isDone(); it.next())
  {
    try {
      #ifdef DEBUG_PRINT
      cout << "Testing: " << endl;
      #endif

      auto gatbIndex = solidKmerOracle->graph->nodeMPHFIndex (it.item());
      #ifdef DEBUG_PRINT
      cout << "gatbIndex = " << gatbIndex << endl;
      cout.flush();
      #endif

      auto gatbKmer = solidKmerOracle->graph->toString (it.item());
      #ifdef DEBUG_PRINT
      cout << "gatbKmer = " << gatbKmer << endl;
      cout.flush();
      #endif

      auto unitigInfo = solidKmerOracle->contains(gatbKmer);
      #ifdef DEBUG_PRINT
      cout << "unitigInfo = " << unitigInfo.toString() << endl;
      cout.flush();
      #endif

      auto unitigSeq = solidKmerOracle->unitigsSequences[unitigInfo.unitigId];
      #ifdef DEBUG_PRINT
      cout << "[before] unitigSeq = " << unitigSeq << endl;
      #endif
      if (unitigInfo.strand=='R')
        unitigSeq = reverse_complement(unitigSeq);
      #ifdef DEBUG_PRINT
      cout << "[after] unitigSeq = " << unitigSeq << endl;
      cout.flush();
      #endif

      if (unitigSeq.find(gatbKmer)!=string::npos) {
        #ifdef DEBUG_PRINT
        cout << "OK!" << endl;
        cout.flush();
        #endif
      }
      else {
        cout << "EYTATester::testSolidKmerOracle() FAILED!!!" << endl;
        cout.flush();
        exit(1);
      }
    }catch (exception& e)
    {
      cout << e.what() << '\n';
      exit(1);
    }
  }

  cout << "EYTATester::testSolidKmerOracle() test - All ok!!!" << endl;
  cout.flush();
}


void EYTATester::checkIfAddingBackTranscriptomicEdgeIsCorrect
                                       (const set<string> &SRGraphRemovedEdges,
                                        const set<string> &edgesInTheTranscriptomicGraph,
                                        const set<string> &correctlyRemovedEdges) {
  cerr << "EYTATester::checkIfEdgeRemovalIsOk() begin test..." << endl;
  cerr.flush();

  for (const string& SRGraphRemovedEdge : SRGraphRemovedEdges) {
    if (edgesInTheTranscriptomicGraph.find(SRGraphRemovedEdge)!=edgesInTheTranscriptomicGraph.end()) //edge is in the transcriptomic graph
    {
      //then it should not be in the correctlyRemovedEdges
      if (correctlyRemovedEdges.find(SRGraphRemovedEdge)!=correctlyRemovedEdges.end()) {
        cerr << "Error in EYTATester::checkIfEdgeRemovalIsOk(): edge " << SRGraphRemovedEdge << " is in edgesInTheTranscriptomicGraph, but also in correctlyRemovedEdges." << endl;
        exit(1);
      }
    }else {
      //edge is not in the transcriptomic graph
      //then it should be in the correctlyRemovedEdges
      if (correctlyRemovedEdges.find(SRGraphRemovedEdge)==correctlyRemovedEdges.end()) {
        cerr << "Error in EYTATester::checkIfEdgeRemovalIsOk(): edge " << SRGraphRemovedEdge << " is NOT in edgesInTheTranscriptomicGraph, and also NOT in correctlyRemovedEdges." << endl;
        exit(1);
      }
    }
  }

  cerr << "EYTATester::checkIfEdgeRemovalIsOk() test - All ok!!!" << endl;
  cerr.flush();
}


void EYTATester::checkRemoveEdgesFromFile(Graph &graph, const string &prefixSR, const string &prefixTr, int nbCores, const set<string> &correctlyRemovedEdges) {
  cerr << "Starting EYTATester::checkRemoveEdgesFromFile..." << endl;
  {
    ofstream debugFile;
    openFileForWriting("debug_removeEdgesFromFile", debugFile);
    graph.removeEdgesFromFile(prefixSR+"_"+prefixTr+".removed_edges", nbCores, &debugFile);
    debugFile.close();
  }

  vector<string> debugVectorString = getVectorStringFromFile("debug_removeEdgesFromFile");
  set<string> debugSetString(debugVectorString.begin(), debugVectorString.end());
  if (correctlyRemovedEdges != debugSetString) {
    cout << "WARNING: FATAL ERROR ON graph->removeEdgesFromSetOfString(): " << __FILE__ << ":" << __LINE__ << endl;
    exit(1);
  }
  cerr << "EYTATester::checkRemoveEdgesFromFile - All ok!" << endl;
  exit(0);
}
