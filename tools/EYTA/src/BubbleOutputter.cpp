//
// Created by Leandro Ishi Soares de Lima on 19/04/17.
//

#include "BubbleOutputter.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

namespace EnhanceTranscriptome {
    string BubbleOutputter::eventTypeAsStr[]={"SNP/1-baseIndel","MSNPs/LargeIndels","Repeats", "AS"};
    string BubbleOutputter::fileType2Description[]={"blat","bubbles.detailed","blat.clustered"};
    string BubbleOutputter::eventType2Description[]={"SNPs_or_1_base_indels","Multiple_SNPs_or_indels","putative_repeats", "alternative_splicing"};



    bool BubbleOutputter::tooSimilar(const string &s, const string &t, double threshold) {
      int ed = computeEditDistance(s,t);
      return (((double)ed)/max(s.size(), t.size()))<threshold; //TODO: should I change this for min, instead of max?
    }


    //this is like BubbleOutputter::tooSimilar, but should be used when one of the strings can be a lot longer than the other one (i.e. we try to fit one string in the other)
    bool BubbleOutputter::tooSimilarSemiGlobal(const string &s, const string &t, double threshold) {
      int ed = computeEditDistanceOfAShortStringInsideALongerString(s,t);
      return (((double)ed)/min(s.size(), t.size()))<threshold;
    }


    //given a string representing a fasta file:
    //>header
    //[sequence]
    //cluster the sequences and return a subset of the original sequences with copies suppressed
    //These sequences are ordered in order of importance (i.e. if the 1st, 3rd and 8th sequence falls on the same cluster, just print the 1st)
    //In my case, the order of importance is given by the quantification
    string BubbleOutputter::clusterSequences (const string &original, double threshold) {
      //parse the original sequence
      vector<pair<string, string> > sequences;
      {
        vector<string> lines;
        boost::split(lines, original, boost::is_any_of("\n"));
        for (int i = 0; i+1 < lines.size(); i += 2)
          sequences.push_back(make_pair(lines[i], lines[i + 1]));
      }


      //do the clustering
      vector<bool> alreadyProcessed(sequences.size(), false);
      stringstream ss;
      for (int i=0;i<sequences.size();i++) {
        if (alreadyProcessed[i])
          continue;

        //output this sequence
        ss << sequences[i].first << endl << sequences[i].second << endl;
        alreadyProcessed[i]=true;

        //find all sequences similar to this one and mark them
        for (int j=i+1;j<sequences.size();j++) {
          //if (tooSimilar(sequences[i].second, sequences[j].second, threshold))
          if (tooSimilarSemiGlobal(sequences[i].second, sequences[j].second, threshold))
            alreadyProcessed[j]=true;
        }
      }

      return ss.str();
    }

}