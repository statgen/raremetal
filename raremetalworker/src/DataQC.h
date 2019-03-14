////////////////////////////////////////////////////////////////////// 
// DatQC.h 
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

#ifndef __DATAQC_H__
#define __DATAQC_H__

#include "StringHash.h"
#include "Pedigree.h"

class SanityCheck
{
public:
    SanityCheck();

    void Check(Pedigree &ped, FILE *log);

    StringIntHash skippedSNPs;
    StringArray chromosomeVCF;
};

#endif


