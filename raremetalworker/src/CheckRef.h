////////////////////////////////////////////////////////////////////// 
// CheckRef.h 
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

#ifndef __CHECKREF_H__
#define __CHECKREF_H__


class CheckRef
{
public:
    CheckRef();

    static String refFile;

    void CheckReferenceGenome(String &chr, Int &position, String &allele1, String &allele2);
};

#endif


