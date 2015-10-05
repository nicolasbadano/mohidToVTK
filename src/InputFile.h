#pragma once

#include <fstream>
#include <string.h>
#include <vector>
#include <sstream>
#include <iostream>


class InputFile
{
private:
    std::ifstream*          ifsPtr;
    int                     lineNumber;
    char                    lineBuffer[5000];
    bool                    openedFile;
    std::string             commentPrefix;
    std::string             fileName_;

public:

    InputFile(const char * fileName)
    {
        fileName_ = std::string(fileName);
        lineNumber = 0;
        openedFile = false;
        commentPrefix = std::string("#");

        ifsPtr = new std::ifstream(fileName);
        openedFile = true;

    }

    InputFile(const std::string fileName)
    {
        fileName_ = fileName;
        lineNumber = 0;
        openedFile = false;
        commentPrefix = std::string("#");

        ifsPtr = new std::ifstream(fileName.data());
        openedFile = true;

    }

    ~InputFile(void)
    {
        if (openedFile) ifsPtr->close();
        delete ifsPtr;
    }

    int returnTagVector(const char *line, std::vector<std::string> &tagVec)
    {
        int inicio;
        int posActual;
        bool posTermino = false;
        int i = 0;
        char tag[200];

        posTermino = true;
        posActual = 0;
        inicio = 0;

        tagVec.resize(0);

        do {
            if (line[i] != '\t' && line[i] != ' ' && line[i] != '\0') {
                //Tenemos un caracter ==> está por comenzar algo??
                if (posTermino == true) {   //Comienza algo
                    posTermino = false;
                    posActual++;
                    inicio = i;
                }
            } else {    //Tenemos un espacio ==> algo terminó?
                if (posTermino == false) {  //Algo terminó
                    posTermino = true;
                    tagVec.resize(posActual);

                    strncpy(tag, line+inicio, i-inicio);
                    tag[i-inicio] = '\0';
                    std::stringstream   tagStream;
                    tagStream.clear();
                    tagStream << tag;
                    tagVec[posActual-1] = tagStream.str();
                }
            }

            i++;
            //std::cout << "Caracter " << i-1 << line[i-1] << "\n";
        } while (line[i-1] != '\0' && line[i-1] != commentPrefix.data()[0]);

        // Se terminó sin devolver la etiqueta --> error!!
        return 0;
    }

    void stripLeadingAndTrailingBlanks(std::string& StringToModify)
    {
        if(StringToModify.empty()) return;

        int startIndex = StringToModify.find_first_not_of(" ");
        if (startIndex == -1) return;
        int endIndex = StringToModify.find_last_not_of(" ");
        if (endIndex == -1) return;
        std::string tempString = StringToModify;
        StringToModify.erase();

        StringToModify = tempString.substr(startIndex, (endIndex-startIndex+ 1) );
    }

    bool returnNextLineOfTags(std::vector<std::string> &tagVec)
    {
        tagVec.clear();

        if (!openedFile) return false;

        while (true) {
            if (!ifsPtr->getline(lineBuffer, 5000, '\n')) {
                return false;
            }

            lineNumber++;
            //PRINTVARIABLE( lineBuffer );

            // Strip white spaces
            std::string lineStr(lineBuffer);
            stripLeadingAndTrailingBlanks(lineStr);
            //PRINTVARIABLE( lineStr );
            // If it is a comment
            if (commentPrefix.compare(0, commentPrefix.length(), lineStr) != 0) {
                returnTagVector(lineBuffer, tagVec);
                //PRINTVARIABLE( tagVec.size() );
                if (tagVec.size() > 0) {
                    //PRINTVARIABLE( tagVec[0] );
                    return true;
                }
            }
        }
    }

    bool returnNextLine(std::string &line)
    {
        line.clear();

        if (!openedFile) return false;

        while (true) {
            if (!ifsPtr->getline(lineBuffer, 5000, '\n'))
                return false;

            lineNumber++;

            // Strip white spaces
            std::string lineStr(lineBuffer);
            stripLeadingAndTrailingBlanks(lineStr);
            // If it is a comment
            int comp = lineStr.compare(0, commentPrefix.length(), commentPrefix);
            if (comp != 0) {
                std::stringstream   lineStream;
                lineStream << lineBuffer;
                line = lineStream.str();
                if (line.length() > 0) {
                    return true;
                } else {
                    continue;
                }
            }
        }
    }

    void close()
    {
        if (openedFile) ifsPtr->close();
        openedFile = false;
    };

    int getLineNumber()
    {
        return lineNumber;
    };

    std::string getLineInfo()
    {
        std::stringstream lineInfoStream;
        lineInfoStream << fileName_ << ":" << lineNumber;
        return lineInfoStream.str();
    };

    void setCommentPrefix(const std::string &prefix)
    {
        commentPrefix = prefix;
    };

    bool isFileOpen()
    {
        return openedFile;
    };
};

