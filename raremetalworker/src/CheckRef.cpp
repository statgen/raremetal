#include "CheckRef.h"
#include "StringBasics.h"
#include "Error.h"
#include "PreMeta.h"

CheckRef::CheckRef() {}

String CheckRef::refFile = "human.g1k.v37-bs.umfa"; 

void ChckRef::CheckReferenceGenome(String & chr, Int & position, String & allele1, String & allele2)
{


}

void CheckRef::Setup()
{
   GenomeSequence reference(refFile);
   String prevChr;
   String refStr;
   genomeIndex_t chrStart = 0;
}
