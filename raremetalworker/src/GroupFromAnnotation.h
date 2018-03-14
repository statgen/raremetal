#ifndef __GROUPFROMANNOTATION_H__
#define __GROUPFROMANNOTATION_H__

#include "IntArray.h"
#include "StringBasics.h"
#include "StringArray.h"
#include "Constant.h"
#include "StringHash.h"
#include "QuickIndex.h"

class GroupFromAnnotation
{
public:
    GroupFromAnnotation();

    ~GroupFromAnnotation();

    //Input/Output options
    static String vcfInput;
    static String groupFile;
    static String function;
    static String mapFile;
    static bool labelHits;

    StringArray *SNPlist;
    StringArray *SNPNoAllele;
    StringArray annoGroups;
    StringArray chrom;
    //these are for annotating single variant accoding to refFlat_hg19.txt
    IntArray start_pos, end_pos;
    StringArray genename, chr;
    StringIntHash ChrStartHash, ChrEndHash;
    QuickIndex chr_idx;

    //this is the position of each marker in a gene for output

    void GetGeneMap(String path);

    void GetGroupFromFile();

    void GetGroupFromVCF();

    void Run(String path);

    String AnnotateSingleVar(String chr, int pos);
};

#endif
