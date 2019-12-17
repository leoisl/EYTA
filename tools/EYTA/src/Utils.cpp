//
// Created by Leandro Ishi Soares de Lima on 05/06/16.
//

#include "Utils.h"
#include <ctime>
#include <sstream>
#include <iostream>
#include <limits>

using namespace std;

char complement(char b)
{
  switch(b)
  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';

    case 'a': return 't';
    case 't': return 'a';
    case 'g': return 'c';
    case 'c': return 'g';

    case 'N': return 'N';
    case '*': return '*';
  }
  return '?';
}

string reverse_complement(const string &seq)
{
  string s(seq.begin(),seq.end());
  string::iterator pos;

  reverse(s.begin(), s.end());

  for(pos=s.begin();pos!=s.end();++pos)
    *pos=complement(*pos);

  return s;
}

//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getAllReadFilesNames(const string &readsFile) {
  vector<string> allReadFilesNames;
  string tempStr;

  ifstream readsFileStream;
  openFileForReading(readsFile, readsFileStream);
  while (getline(readsFileStream, tempStr)) {
    if (tempStr.size() > 0)
      allReadFilesNames.push_back(tempStr);
  }
  readsFileStream.close();

  return allReadFilesNames;
}

void checkIfFileExists(const string &filename) {
  ifstream input;
  openFileForReading(filename, input);
  input.close();
}

//check if a read file exists and all files inside it
//dies to a fatal error if any of the files do not exist
void checkIfReadFileIsFine(const string &readsFileName) {
  vector<string> allReadsFiles = getAllReadFilesNames(readsFileName);
  for (const string &readFileName : allReadsFiles)
    checkIfFileExists(readFileName);
}

int computeEditDistanceCore(const string &s1, const string &s2) {
  const std::size_t len1 = s1.size(), len2 = s2.size();
  std::vector<int> col(len2+1), prevCol(len2+1);

  for (int i = 0; i < prevCol.size(); i++)
    prevCol[i] = i;
  for (int i = 0; i < len1; i++) {
    col[0] = i+1;
    for (int j = 0; j < len2; j++)
      // note that std::min({arg1, arg2, arg3}) works only in C++11,
      // for C++98 use std::min(std::min(arg1, arg2), arg3)
      //TODO: treat the evetuel Ns that we can have here?
      col[j+1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i]==s2[j] ? 0 : 1) });
    col.swap(prevCol);
  }
  return prevCol[len2];
}

int computeEditDistance(const string &s1, const string &s2) {
  return min(computeEditDistanceCore(s1, s2), computeEditDistanceCore(reverse_complement(s1), s2));
}



//This is a function to compute the minimum edit distance of a short string inside a longer string
//YOU SHOULD NOT CALL THIS FUNCTION, BUT computeEditDistanceOfAShortStringInsideALongerString()
//Basically, the idea is:
//You have:
//A long string:  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//A short string: yyyyyyyyyyyyyyyyyyyy
//And you want to know the best edit distance of the short string to a substring of the long string
//In order to do so, we compute the edit distance in a common way, where the short read is in x-axis, long read is in y-axis,
//and the first row is all 0s
//and the result is the best edit distance in the last row
int computeEditDistanceOfAShortStringInsideALongerStringCore(const string &longString, const string &shortString) {
  const std::size_t len1 = longString.size(), len2 = shortString.size();
  std::vector<int> col(len2+1), prevCol(len2+1);

  //fills the first collumn
  for (int i = 0; i < prevCol.size(); i++)
    prevCol[i] = i;

  int editDistance = prevCol.back();
  int bestIndex=0;

  for (int i = 0; i < len1; i++) {
    col[0] = 0; //first row is all 0s
    for (int j = 0; j < len2; j++)
      // note that std::min({arg1, arg2, arg3}) works only in C++11,
      // for C++98 use std::min(std::min(arg1, arg2), arg3)
      //TODO: treat the evetuel Ns that we can have here?
      col[j+1] = std::min({ prevCol[j+1] + 1, col[j] + 1, prevCol[j] + (longString[i]==shortString[j] ? 0 : 1) });

    if (col.back()<editDistance) {
      editDistance = col.back();
      bestIndex = i;
    }

    //editDistance=std::min(editDistance, col.back());
    col.swap(prevCol);
  }

  return editDistance;
}


int computeEditDistanceOfAShortStringInsideALongerString(const string &s1, const string &s2) {
  string longString, shortString;
  if (s1.size() < s2.size()) {
    shortString = s1;
    longString = s2;
  }else {
    shortString = s2;
    longString = s1;
  }

  return min(computeEditDistanceOfAShortStringInsideALongerStringCore(longString, shortString),
             computeEditDistanceOfAShortStringInsideALongerStringCore(reverse_complement(longString), shortString));
}



//find the best anchors of the kmer into the longString
vector<size_t> getBestAnchors(const string &longString, const string &kMer) {
  vector<size_t> anchors;
  int minEditDistance=numeric_limits<int>::max();

  for (size_t i=0;i<longString.size()-kMer.size(); i++) {
    auto lsKmer = longString.substr(i, kMer.size());
    int editDistance=computeEditDistanceCore(lsKmer, kMer);
    if (editDistance < minEditDistance) {
      minEditDistance = editDistance;
      anchors.clear();
      anchors.push_back(i);
    }else if (editDistance == minEditDistance) {
      anchors.push_back(i);
    }
  }

  return anchors;
}

//compute the "anchored distance" of the shortest string to the longest string
//The idea is that the k-1 first characters and the k-1 last characters should be in the longest string
//Thus, we find the positions where the k-1 first characters best anchor in the longest string (i.e. for each (k-1)-mer of the longest string, we compute the edit distance of the (k-1)-mer to the first (k-1)-mer of the short string. We record the best positions)
//Do the same for the last (k-1)-mer of the shortest string
//For each pair of positions i, j, where i<j and i is the start position of a best anchor of the first (k-1)-mer and j is the end position of a best anchor of the last (k-1)-mer
    //get the LR string between i and j, and compute the ED between this substring and the short string
    //The anchored distance is the least distance of the ED distance computed in the previous step
int computeAnchoredDistanceCore(const string &longString, const string &shortString, int k) {
  k=k-1; //easier to work

  //compute the anchors
  auto firstKmerSS = shortString.substr(0, k);
  auto firstKmerSSAnchors = getBestAnchors(longString, firstKmerSS);
  auto lastKmerSS = shortString.substr(shortString.size()-k, k);
  auto lastKmerSSAnchors = getBestAnchors(longString, lastKmerSS);

  //compute the minimum anchored distance
  int minAnchoredDistance=numeric_limits<int>::max();
  for(const auto &firstKmerAnchor : firstKmerSSAnchors) {
    for(const auto &lastKmerAnchor : lastKmerSSAnchors) {
      int editDistance=computeEditDistanceCore(longString.substr(firstKmerAnchor, lastKmerAnchor+k-firstKmerAnchor), shortString);
      if (editDistance<minAnchoredDistance)
        minAnchoredDistance=editDistance;
    }
  }

  return minAnchoredDistance;
}

int computeAnchoredDistance(const string &longString, const string &shortString, int k) {
  return min(computeAnchoredDistanceCore(longString, shortString, k),
             computeAnchoredDistanceCore(reverse_complement(longString), shortString, k));
}


string getTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[1024];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}



string toUpperCase (const string &s) {
  string toReturn(s);
  //transform the string
  for (int j=0;j<s.size();j++)
    toReturn[j]=toupper(s[j]);
  return toReturn;
}


//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile) {
  vector<string> allReadFilesNames;
  string tempStr;

  ifstream readsFileStream;
  openFileForReading(readsFile, readsFileStream);
  while (getline(readsFileStream, tempStr)) {
    if (tempStr.size() > 0)
      allReadFilesNames.push_back(tempStr);
  }
  readsFileStream.close();

  return allReadFilesNames;
}

void openFileForReading(const string &filePath, ifstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void openFileForWriting(const string &filePath, ofstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void fatalError (const string &message) {
  cerr << endl << endl << "[FATAL ERROR] " << message << endl << endl;
  cerr.flush();
  exit(1);
}


long countNbNodes (const string &unitigSizeFilename) {
  ifstream unitigSizeFile;
  openFileForReading(unitigSizeFilename, unitigSizeFile);
  long nbNodes;
  unitigSizeFile >> nbNodes;
  unitigSizeFile.close();
  return nbNodes;
}


//read removed edges into a set of string
set<string> readRemovedEdgesIntoASetOfString(const string &filename) {
  set<string> removedEdges;
  vector<string> removedEdgesAsVector = getVectorStringFromFile(filename);
  removedEdges.insert(removedEdgesAsVector.begin(), removedEdgesAsVector.end());
  return removedEdges;
}



void concatenate2FilesIntoA3rd(const string &filename1, const string &filename2, const string &filename3) {
  ifstream file1;
  ifstream file2;
  ofstream file3;

  openFileForReading(filename1, file1);
  openFileForReading(filename2, file2);
  openFileForWriting(filename3, file3);

  //concat
  string line;
  while(getline(file1, line))
    file3 << line << endl;
  while(getline(file2, line))
    file3 << line << endl;

  file1.close();
  file2.close();
  file3.close();
}