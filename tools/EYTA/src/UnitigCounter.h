//
// Created by Leandro Ishi Soares de Lima on 15/11/17.
//

#ifndef KSGATB_UNITIGCOUNTER_H
#define KSGATB_UNITIGCOUNTER_H

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/
#include <set>

using namespace std;

//Represents an unitig counter associated to a color (blue if transcriptome, red if short reads):
// 1. Counts how many reads map to a given unitig
class UnitigCounter {
private:
    vector<long> count; //the counts - the index is the unitig index
    char color; //the color of a node if the counts of that node is > 0
    ISynchronizer* synchro; //controls multithreaded access to addToCount() function
public:
    //constructor
    UnitigCounter (long nbUnitigs, char color) : count(nbUnitigs), color(color), synchro(System::thread().newSynchronizer()){}

    //add count to a node
    void addToCount(long node, long amount) {
      synchro->lock();
      count[node]+=amount;
      synchro->unlock();
    }

    //add count to several nodes
    //more efficient than calling addToCount() several times in multi-threaded environments
    void addToCount(const set<long> &nodes, long amount) {
      synchro->lock();
      for (long node : nodes)
        count[node]+=amount;
      synchro->unlock();
    }

    //get the node color
    string getNodeColor(long node) const {
      string toReturn="";
      if (count[node]>0)
        toReturn+=color;
      return toReturn;
    }

    //get the node count
    long getNodeCount(long node) const {
      return count[node];
    }

    long getInvalidUnitigId() const {
      return -1;
    }
};


#endif //KSGATB_UNITIGCOUNTER_H
