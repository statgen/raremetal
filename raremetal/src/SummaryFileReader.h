////////////////////////////////////////////////////////////////////// 
// SummaryFileReader.h 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 

#ifndef __SUMMARYFILEREADER_H__
#define __SUMMARYFILEREADER_H__

#include "Tabix.h"
#include "InputFile.h"
#include "StringBasics.h"
#include "StringArray.h"
#include "StringHash.h"

class SummaryFileReader
{
public:
    static int counter;
    int marker_col;
    int cov_col;

    SummaryFileReader()
    {
        old_chr = "";
        old_pos = 0;
        counter = 0;
    };

    String marker_cov, marker_nearby, old_chr;
    StringArray markerNearby, markerNearbyCov;
    StringIntHash markerPosHash;
    String buffer;
    int old_pos;

    Tabix myTabix;
    IFILE myFilePtr;

    bool ReadTabix(String filename);

    bool ReadRecord(String chr, int pos);
};

#endif


