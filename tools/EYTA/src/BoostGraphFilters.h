//
// Created by Leandro Ishi Soares de Lima on 18/04/17.
//

#ifndef KSGATB_BOOSTGRAPHFILTERS_H
#define KSGATB_BOOSTGRAPHFILTERS_H

#include <vector>
#include <set>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "EnhanceTranscriptomeDefs.h"

using namespace std;

namespace EnhanceTranscriptome {
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //class related to the vertex and edge filters used when defining the subgraphs
    //This class uses a vector<bool> as predicate: i.e. its size is the number of nodes, but the queries as O(1)
    class SubgraphEdgeFilter {
    private:
        graph_t *g;
        vector<bool> predicate; //TODO: Memory issue: change this for boost::dynamic_bitset or something like this
    public:
        SubgraphEdgeFilter() : g(NULL), predicate() { } //default constructor, needed by boost
        SubgraphEdgeFilter(graph_t *g) : g(g), predicate(num_edges(*g), true) { } //at the start, the graph is full!

        //TODO: using a validating cycle method (see cormen-solutions, page 160), we don't need reset function - all operations will be O(1)
        //TODO: this function is called a lot, so I think it is valid to do this...
        void reset(bool value = false) {
          for (size_t i = 0; i < predicate.size(); ++i)
            predicate[i] = value;
        }
        //TODO: this function is called a lot, so I think it is valid to do this...
        //TODO: using a validating cycle method (see cormen-solutions, page 160), we don't need reset function - all operations will be O(1)

        void add(const Edge &e) {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return; //artificial edges are not indexed here
          predicate[id] = true;
        }

        void remove(const Edge &e) {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return; //artificial edges are not indexed here
          predicate[id] = false;
        }

        //check if the edge is present or not
        bool check(const Edge &e) const {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return true; //artificial edges are always present
          return predicate[id];
        }
    };

    //This class uses a vector<bool> as predicate: i.e. its size is the number of nodes, but the queries as O(1)
    class SubgraphVertexFilter {
    private:
        graph_t *g;
        vector<bool> predicate; //TODO: Memory issue: change this for boost::dynamic_bitset or something like this
        SubgraphEdgeFilter *subgraphEdgeFilter;
    public:
        SubgraphVertexFilter() : g(NULL), predicate(),
                                 subgraphEdgeFilter(NULL) { } //default constructor, needed by boost
        SubgraphVertexFilter(graph_t *g, SubgraphEdgeFilter *subgraphEdgeFilter) : g(g),
                                                                                   predicate(num_vertices(*g), true),
                                                                                   subgraphEdgeFilter(
                                                                                       subgraphEdgeFilter) { } //at the start, the graph is full!

        //TODO: using a validating cycle method (see cormen-solutions, page 160), we don't need reset function - all operations will be O(1)
        //TODO: this function is called a lot, so I think it is valid to do this...
        void reset(bool value = false) {
          subgraphEdgeFilter->reset(value);
          for (size_t i = 0; i < predicate.size(); ++i)
            predicate[i] = value;
        }
        //TODO: this function is called a lot, so I think it is valid to do this...
        //TODO: using a validating cycle method (see cormen-solutions, page 160), we don't need reset function - all operations will be O(1)

        void add(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          predicate[boost::get(boost::vertex_index, *g, v)] = true;
        }

        void remove(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          predicate[boost::get(boost::vertex_index, *g, v)] = false;
        }

        //check if the vertex is present or not
        bool check(const Vertex &v) const {
          return predicate[boost::get(boost::vertex_index, *g, v)];
        }


        void removeAllBlueOnlyNodes() {
          for (size_t i = 0; i < predicate.size(); ++i) {
            Vertex v = boost::vertex(i, *g);
            if ((*g)[v].transcriptomeCoverage > 0 && (*g)[v].shortReadsCoverage==0) //the node is blue-only
              remove(v); //remove it
          }
        }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////









    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //This class uses a set<Edge> as predicate, where if the edge is in the set, then the edge is removed
    //Good to use when you remove VERY FEW things from the graph
    class SubgraphEdgeFilterRemoveSet {
    private:
        graph_t *g;
        set<Edge> predicate;
    public:
        SubgraphEdgeFilterRemoveSet() : g(NULL), predicate() { } //default constructor, needed by boost
        SubgraphEdgeFilterRemoveSet(graph_t *g) : g(g), predicate() { } //at the start, the graph is full!

        void putEveryoneBack() { predicate.clear(); }

        void add(const Edge &e) {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return; //artificial edges are not indexed here
          predicate.erase(e);
        }

        void remove(const Edge &e) {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return; //artificial edges are not indexed here
          predicate.insert(e);
        }

        //check if the edge is present or not
        bool check(const Edge &e) const {
          int id = (*g)[e].id;
          if (id == ARTIFICIAL_ID) return true; //artificial edges are always present
          return predicate.find(e)==predicate.end();
        }
    };

    //This class uses a set<Vertex> as predicate, where if the edge is in the set, then the edge is removed
    //Good to use when you remove VERY FEW things from the graph
    class SubgraphVertexFilterRemoveSet {
    private:
        graph_t *g;
        set<Vertex> predicate;
        SubgraphEdgeFilterRemoveSet *subgraphEdgeFilter;
    public:
        SubgraphVertexFilterRemoveSet() : g(NULL), predicate(),
                                 subgraphEdgeFilter(NULL) { } //default constructor, needed by boost
        SubgraphVertexFilterRemoveSet(graph_t *g, SubgraphEdgeFilterRemoveSet *subgraphEdgeFilter) : g(g),
                                                                                   predicate(),
                                                                                   subgraphEdgeFilter(
                                                                                       subgraphEdgeFilter) { } //at the start, the graph is full!

        void putEveryoneBack() {
          subgraphEdgeFilter->putEveryoneBack();
          predicate.clear();
        }

        void add(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          predicate.erase(v);
        }

        void remove(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          predicate.insert(v);
        }

        //check if the vertex is present or not
        bool check(const Vertex &v) const {
          return predicate.find(v)==predicate.end();
        }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////



    //This is what has to be given to the Filtered graph
    template <class VertexFilterClass>
    class SubgraphVertexFilterForFG {
    public:
        VertexFilterClass *subgraphVertexFilter; //An object simply does not work here... It needs to be a pointer... I have no goddamn idea why
        SubgraphVertexFilterForFG():subgraphVertexFilter(NULL){}
        SubgraphVertexFilterForFG(VertexFilterClass *subgraphVertexFilter):subgraphVertexFilter(subgraphVertexFilter){}
        bool operator()(const Vertex &v) const { return subgraphVertexFilter->check(v); }
    };
    //This is what has to be given to the Filtered graph
    template <class EdgeFilterClass>
    class SubgraphEdgeFilterForFG {
    public:
        EdgeFilterClass *subgraphEdgeFilter; //An object simply does not work here... It needs to be a pointer... I have no goddamn idea why
        SubgraphEdgeFilterForFG():subgraphEdgeFilter(NULL){}
        SubgraphEdgeFilterForFG(EdgeFilterClass *subgraphEdgeFilter):subgraphEdgeFilter(subgraphEdgeFilter){}
        bool operator()(const Edge &e) const { return subgraphEdgeFilter->check(e); }
    };
    typedef boost::filtered_graph <graph_t, SubgraphEdgeFilterForFG<SubgraphEdgeFilter>, SubgraphVertexFilterForFG<SubgraphVertexFilter> > FilteredGraph;
    typedef boost::graph_traits<FilteredGraph>::adjacency_iterator AdjacencyIterator;
    typedef boost::filtered_graph <graph_t, SubgraphEdgeFilterForFG<SubgraphEdgeFilterRemoveSet>,
        SubgraphVertexFilterForFG<SubgraphVertexFilterRemoveSet> > FilteredGraphRemoveSet;
}

#endif //KSGATB_BOOSTGRAPHFILTERS_H
