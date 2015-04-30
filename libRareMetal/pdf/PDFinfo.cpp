////////////////////////////////////////////////////////////////////// 
// janpdf/PDFinfo.cpp 
// (c) 2000-2007 Goncalo Abecasis (c) 2002-2007 Jan Wigginton
// 
// This file is distributed as part of the PEDSTATS source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile PEDSTATS.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#include "PDFinfo.h"
#include "PDF.h"

PDFInfo::PDFInfo()
   {
   time_t creation = time(NULL);
   tm * time_struct = gmtime(&creation);

   creationYear = time_struct->tm_year + 1900;
   creationMonth = time_struct->tm_mon + 1;
   creationDay = time_struct->tm_mday;

   creator = "A C++ PDF Library (c) 2001 Goncalo Abecasis";
   }

void PDFInfo::Write(PDF & pdf)
   {
   object = pdf.GetObject();

   pdf.OpenObject(object);
   pdf.OpenDictionary();

   if (!title.IsEmpty()) pdf.WriteString("Title", title);
   if (!author.IsEmpty()) pdf.WriteString("Author", author);
   if (!subject.IsEmpty()) pdf.WriteString("Subject", subject);
   if (!keywords.IsEmpty()) pdf.WriteString("Keywords", keywords);
   if (!producer.IsEmpty()) pdf.WriteString("Producer", producer);

   pdf.WriteString("Creator", creator);
   pdf.WriteDate("CreationDate", creationYear, creationMonth, creationDay);

   pdf.CloseDictionary();
   pdf.CloseObject();
   } 
 
 
 
 
 
 
 
 
 
