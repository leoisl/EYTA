//
// Created by Leandro Ishi Soares de Lima on 18/04/17.
//

#ifndef KSGATB_ENHANCETRANSCRIPTOMEDEFS_H
#define KSGATB_ENHANCETRANSCRIPTOMEDEFS_H

#define ARTIFICIAL_ID -322 //any negative value in fact

#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/scope_exit.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace std;

namespace EnhanceTranscriptome {
    //some definitions used by everyone in EnhanceTranscriptome Module

    //define the structures needed by these classes
    //vertex informations
    struct VertexInfo {
        int id; //do not use unsigned values
        bool isTrustable; //is this unitig trustable or not?
        int transcriptomeCoverage;
        int shortReadsCoverage;
        string seq; //probably represent this as (unitig id, pos)
        char strand; //strand
        int weight; //how many kmers we have in this unitig?
    };

    //edge informations
    struct EdgeInfo {
        int id; //do not use unsigned values
        int sourceWeight; //the weight of the source unitig
        int targetWeight; //the weight of the target unitig
        int edgeWeight; //how many kmers does this edge represent?
        string seq; //the sequence of this edge
    };

    //some typedefs for making life easier
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::bidirectionalS, VertexInfo, EdgeInfo> graph_t;
    typedef boost::graph_traits<graph_t>::vertex_descriptor Vertex;
    typedef boost::graph_traits<graph_t>::edge_descriptor Edge;

    //Type used to represent a path.
    //It keeps track of the length (distance) of the path, its number of branching nodes and its nodes
    template <class GraphType>
    class Path {
    public:
        //TODO: move this to private with proper getters and setters
        int distance;
        int branchingNodes;
        vector<Vertex> nodes;
        const GraphType *graph; //A path belongs to a Graph. this denotes the graph this path belongs to.

        Path(const GraphType *graph):distance(0), branchingNodes(0), nodes(), graph(graph){}
        Path (int distance, int branchingNodes, const vector<Vertex> &nodes, const GraphType *graph) : distance(distance), branchingNodes(branchingNodes), nodes(nodes), graph(graph){}

        double getQuantification() const {
          if (nodes.size()==0) return 0;

          double quant=0;
          for (const auto &node : nodes)
            quant += (*graph)[node].shortReadsCoverage;
          return quant/nodes.size();
        }

        string getNodesAsStr () const {
          stringstream ss;
          for(const auto &node : nodes)
            ss << (*graph)[node].id << "_" << (*graph)[node].strand << "_";
          return ss.str();
        }
    };


    //mapping info
    class MappingInfo {
    public:
        static vector<MappingInfo> readAllMappingInfo (const string &mappingReadsFilename, graph_t& graph, const map<pair<int, char>, int> &labelToIndex);
        static void readTranscriptSequences (const string &mappingReadsFilename, map<long, string> &transcriptIndex2sequence);
        const vector<Vertex> & getMappingFW() const { return mappingFW; }

        //TODO: do we need this???
        const vector<Vertex> & getMappingRC() const { return mappingRC; }
        //TODO: do we need this???

        const set<Vertex> & getUniqueNodes() const { return uniqueNodes; }
        const string & getRead() const { return read; }
        const string & getHeader() const { return header; }
        double getSRQuantification() const { return SRQuantification; }
        int getReadIndex() const { return readIndex; }
    private:
        MappingInfo(){} //no one can instatiate this, only with readAllMappingInfo()
        int readIndex;
        string header;
        string read;

        //these represent the mappings in the FW and RC
        vector<Vertex> mappingFW;
        vector<Vertex> mappingRC;
        vector<int> positionInsideUnitigFW;
        vector<int> positionInsideUnitigRC;

        //these are the unique nodes a read map to
        set<Vertex> uniqueNodes;

        //quantification of the LRs according to the SRs it maps to
        double SRQuantification;

        //read the mapping
        static void readTheMapping(ifstream& readMappingReader, vector <Vertex>& mapping, vector<int> &positionInsideUnitig, graph_t& graph, const map<pair<int, char>, int> &labelToIndex);
    };

    class GraphWriter {
    public:
        //some prints to help
        template<class Graph>
        static string toString(const Vertex &v, const Graph &graph) {
          stringstream ss;
          ss << graph[v].id << "_" << graph[v].strand;
          return ss.str();
        }

        //some prints to help
        template<class Graph>
        static string toStringNodeComplete(const Vertex &v, const Graph &graph) {
          stringstream ss;
          string color="";
          if (graph[v].transcriptomeCoverage>0) color+='B';
          if (graph[v].shortReadsCoverage>0) color+='R';
          ss << graph[v].id << "_" << graph[v].strand << "\t" << color << "\t" <<  graph[v].seq << "\t" <<  graph[v].weight;
          return ss.str();
        }

        template<class Graph>
        static string toString(const Edge &e, const Graph &graph) {
          stringstream ss;
          ss << toString(source(e, graph), graph) << "\t" << toString(target(e, graph), graph) << "\t" << graph[e].targetWeight;
          return ss.str();
        }

        template <class GraphType>
        static string toString(const Path<GraphType> &path, const GraphType &g) {
          stringstream ss;
          ss << "Path:" << endl;
          ss << "Nodes = ";
          for (const auto& v : path.nodes)
            ss << GraphWriter::toString(v, g) << " ";
          ss << endl << "Distance = " << path.distance << endl << "Branching nodes = " << path.branchingNodes;
          return ss.str();
        }

        //declare the graphviz writer
        //TODO: fix this - it is not working
        //    class graphviz_writer {
        //    public:
        //        graphviz_writer(const graph_t &graph) : graph(graph) { }
        //
        //        template<class AnythingElse>
        //        void operator()(std::ostream &out, const AnythingElse &e) const {
        //          //do nothing for other types
        //        }
        //
        //    private:
        //        const graph_t &graph;
        //    };
        //    //template specialization for Vertex and Edge
        //    template<>
        //    inline void graphviz_writer::operator()(std::ostream &out, const Vertex &v) const {
        //      out << "[label=\"" << graph[v].id << "_" << graph[v].strand << "\"]";
        //    }
        //    template<>
        //    inline void graphviz_writer::operator()(std::ostream &out, const Edge &e) const {
        //      out << "[label=\"" << graph[e].targetWeight << "\"]";
        //    }
        //    //function to write the graph to the file
        //    void writeGraphToFile(const graph_t &graph, const string &fileName) {
        //      //output the graph
        //      ofstream graphOut(fileName);
        //      write_graphviz(graphOut, graph, graphviz_writer(graph));
        //      graphOut.close();
        //    }

        template<class T>
        static void writeGraphToFile(const T &graph, const string &fileName) {
          //output the graph
          ofstream graphOut(fileName);

          typedef typename boost::graph_traits<T>::vertex_iterator vertex_iter;
          std::pair <vertex_iter, vertex_iter> vp;
          for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
            Vertex v = *(vp.first);
            auto outEdges = out_edges(v, graph);
            for_each(outEdges.first, outEdges.second, [&](const Edge &e) {
                graphOut << toString(e, graph) << endl;
            });
          }
          graphOut.close();
        }

        template <class GraphType>
        static string toString(const vector<Vertex> &alternativePath, GraphType &graph) {
          stringstream ss;
          for_each(alternativePath.begin(), alternativePath.end(), [&](const Vertex &v) {
              ss << GraphWriter::toString(v, graph) << " ";
          });
          return ss.str();
        }


        template <class GraphType>
        static string toString(const vector<Vertex> &firstPart, int numberOfNs, const vector<Vertex> &secondPart, GraphType &graph) {
          stringstream ss;
          for_each(firstPart.begin(), firstPart.end(), [&](const Vertex &v) {
              ss << GraphWriter::toString(v, graph) << " ";
          });
          ss << numberOfNs << "Ns ";
          for_each(secondPart.begin(), secondPart.end(), [&](const Vertex &v) {
              ss << GraphWriter::toString(v, graph) << " ";
          });
          return ss.str();
        }

    };



    //classes related to the Dijkstra visitor for finding alternative paths
    struct LongReadWasReached {
        int distance; //distance to the first node in the LR that was reached
        LongReadWasReached(int distance) : distance(distance) { }
    };

    struct LongReadWasNotReached {
    };

    class StopWhenReachingLRDijkstraVisitor : public boost::dijkstra_visitor<> {
    private:
        const set <Vertex> &targetNodes;
        vector<int> &distances;
        int maxDistance;
    public:
        StopWhenReachingLRDijkstraVisitor(const set <Vertex> &targetNodes, vector<int> &distances, int maxDistance) :
            targetNodes(targetNodes), distances(distances), maxDistance(maxDistance) { }

        template<class Vertex, class Graph>
        void finish_vertex(Vertex v, const Graph &g) {
          //cout << "[StopWhenReachingLRDijkstraVisitor] Vertex " << g[v].id << "_" << g[v].strand << " finished with distance " << distances[boost::get(boost::vertex_index, g, v)] << endl;

          //check if we are already going over the max distance - since dijkstra is a greedy algorithm, then we can stop here!
          if (distances[boost::get(boost::vertex_index, g, v)] > maxDistance)
            throw LongReadWasNotReached();

          //check if we reached the target node within the allowed distance
          if (targetNodes.find(v) != targetNodes.end())
            throw LongReadWasReached(distances[boost::get(boost::vertex_index, g, v)]);
        }
    };

    //class related to get the subgraph
    struct TooDistant {
    };

    class StopWhenVeryDistantFromSourceDijkstraVisitor : public boost::dijkstra_visitor<> {
    private:
        vector<int> &distances;
        int maxDistance;
        set <Vertex> &verticesReachableByTheSourceWithinMaxDistance;
    public:
        StopWhenVeryDistantFromSourceDijkstraVisitor(vector<int> &distances, int maxDistance,
                                                     set <Vertex> &verticesReachableByTheSourceWithinMaxDistance) :
            distances(distances), maxDistance(maxDistance),
            verticesReachableByTheSourceWithinMaxDistance(verticesReachableByTheSourceWithinMaxDistance) { }

        template<class Vertex, class Graph>
        void finish_vertex(Vertex v, const Graph &g) {
          //EYTA_DEBUG
          //cout << "[StopWhenVeryDistantFromSourceDijkstraVisitor] Vertex " << toString(v,g) << " finished - distance = " << distances[boost::get(boost::vertex_index, g, v)] << endl;
          //EYTA_DEBUG

          //check if we are already going over the max distance - since dijkstra is a greedy algorithm, then we can stop here!
          if (distances[boost::get(boost::vertex_index, g, v)] > maxDistance)
            throw TooDistant();
          else //this node is reachable!
            verticesReachableByTheSourceWithinMaxDistance.insert(v);
        }
    };

    //helper function that says if a node is branching
    template <class GraphType>
    int isBranching(const Vertex &v, GraphType &g) {
      return boost::out_degree(v,g)>1 || boost::in_degree(v,g)>1;
    }
}


#endif //KSGATB_ENHANCETRANSCRIPTOMEDEFS_H
