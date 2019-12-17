//
// Created by Leandro Ishi Soares de Lima on 05/06/16.
//

#ifndef KISSPLICE_UTILS_H
#define KISSPLICE_UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>

using namespace std;

char complement(char b);
string reverse_complement(const string &seq);

//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getAllReadFilesNames(const string &readsFile);

//checks if a read exists or not
void checkIfFileExists(const string &filename);

//check if a read file exists and all files inside it
//dies to a fatal error if any of the files do not exist
void checkIfReadFileIsFine(const string &readsFileName);

int computeEditDistance(const string &originalSequence, const string &alternativeSequence);
int computeEditDistanceOfAShortStringInsideALongerString(const string &longString, const string &shortString);
int computeAnchoredDistance(const string &longString, const string &shortString, int k);

string toUpperCase (const string &s);

string getTime();

//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile);

void openFileForReading(const string &filePath, ifstream &stream);
void openFileForWriting(const string &filePath, ofstream &stream);
void concatenate2FilesIntoA3rd(const string &filename1, const string &filename2, const string &filename3);

void fatalError (const string &message);

long countNbNodes (const string &unitigSizeFilename);

//read removed edges into a set of string
set<string> readRemovedEdgesIntoASetOfString(const string &filename);


#endif //KISSPLICE_UTILS_H
