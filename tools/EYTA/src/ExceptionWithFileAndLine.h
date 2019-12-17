//
// Created by Leandro Ishi Soares de Lima on 12/10/17.
//

#ifndef KSGATB_EXCEPTIONWITHFILEANDLINE_H
#define KSGATB_EXCEPTIONWITHFILEANDLINE_H

#include <exception>
#include <string>

using namespace std;

class ExceptionWithFileAndLine : public exception{
private:
    string reason;
public:
    ExceptionWithFileAndLine(const string &file, int line) {
      reason = string("Exception at ") + file + string(":") + to_string(line);
    }
    const char* what() const noexcept {
      return reason.c_str();
    }
};


#endif //KSGATB_EXCEPTIONWITHFILEANDLINE_H
