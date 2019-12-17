//
// Created by Leandro Ishi Soares de Lima on 11/01/18.
//

#ifndef EYTA_MAPPINGOUTPUTTER_H
#define EYTA_MAPPINGOUTPUTTER_H

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/
#include <fstream>
#include <string>

using namespace std;

//Just a class to output the results of a mapping in a synchronized way
class MappingOutputter {
private:
    ofstream out;
    ISynchronizer* synchro; //controls multithreaded access to output() function
public:
    //constructor: construct a mapping output that will output to filename
    MappingOutputter(const string &filename) : out(filename), synchro(System::thread().newSynchronizer()){}

    //destructor: close the file
    ~MappingOutputter(){ out.close(); }

    //output a str multithreadely
    void output(const string &str) {
      synchro->lock();
        out << str;
      synchro->unlock();
    }

};


#endif //EYTA_MAPPINGOUTPUTTER_H
