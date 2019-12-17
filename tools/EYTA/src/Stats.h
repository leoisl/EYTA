//
// Created by Leandro Ishi Soares de Lima on 01/06/18.
//

#ifndef EYTA_STATS_H
#define EYTA_STATS_H
#include <map>
#include <gatb/gatb_core.hpp>
#include <string>
#include "Utils.h"


using namespace std;



class Stats {
private:
    map<string, ISynchronizer*> synchros; //controls multi-threaded access to the stats
    map<string, long> stats; //the stats themselves

    void checkIfMetricExists(const string &stat) const {
      //TODO: do this better with a method - I had no internet access and did not remember the map definition
      for (const auto& pair : stats) {
        if (pair.first == stat)
          return;
      }

      throw runtime_error("Stat not found: "+stat);
    }
public:
    Stats(){
      //vector with the stats
      vector<string> statsMetrics={
          "nbOfKmersInReadsThatMappedToTheGraph",
          "nbOfKmersInReadsThatDidNotMapToTheGraph",
          "nbOfKmersInReadsThatMappedToTrustfulUnitigs",
          "nbOfKmersInReadsThatMappedToUntrustfulUnitigs",
          "nbOfULGEdgesNotAddedDueToStartUnitigNotGoingUntilTheEnd",
          "nbOfULGEdgesNotAddedDueToEndUnitigNotStartingInTheBegin",
          "nbOfHintsNotAddedDueToEndUnitigNotStartingInTheBegin",
          "nbOfULGEdgesNotAddedDueToMappingToEqualUnitigs",
          "nbOfULGEdgesAdded",
          "nbOfHintsAdded"
      };

      //stats init
      for (const string &metric : statsMetrics) {
        synchros[metric]=System::thread().newSynchronizer();
        stats[metric]=0;
      }
    }

    void increaseStats(const string &metric, int increaseAmount=1) {
      checkIfMetricExists(metric);

      synchros.at(metric)->lock();
        stats[metric]+=increaseAmount;
      synchros.at(metric)->unlock();
    }

    long getStats(const string &metric) const {
      checkIfMetricExists(metric);

      synchros.at(metric)->lock();
        return stats.at(metric);
      synchros.at(metric)->unlock();
    }

    void printStatsToFile(const string &filename) const {
      ofstream file;
      openFileForWriting(filename, file);
      for (const auto& pair : stats)
        file << pair.first << " = " << getStats(pair.first) << endl;
      file.close();
    }
};


#endif //EYTA_STATS_H
