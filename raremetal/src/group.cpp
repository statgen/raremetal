#include "Parameters.h"
#include "math.h"
#include <iostream>
#include <stdio.h>
#include <time.h>
#include "Error.h"
#include "StringBasic.h"

int main(int argc, char **argv)
{
    String infile, outfile, function;

    BEGIN_LONG_PARAMETERS(additional)
                    LONG_PARAMETER_GROUP("Input/output:")
                    LONG_STRINGPARAMETER("in", &infile)
                    LONG_STRINGPARAMETER("out", &outfile)
                    LONG_STRINGPARAMETER("function", &function)

    END_LONG_PARAMETERS();


    ParameterList pl;
    pl.Add(new LongParameters("Options:", additional));

    pl.Read(argc, argv);
    pl.Status();

    FILE *in;
    in = fopen(infile, "r");
    String buffer;
    StringArray tmp;
    int genecount = 0;

    while (!feof(inFile))
    {
        buffer.ReadLine(inFile);
        tmp.Clear();
        tmp.AddTokens(buffer, "\t;");
        for (int i = 0; i < tmp.Length(); i++)
        {
            if (tmp[i].Find(function) != -1)
            {
                StringArray sub;
                sub.AddTokens(tmp[i], ":");
                String genename = sub[1];
                int gene = groupHash.Integer(genename);
                if (gene == -1)
                {
                    groupHash.SetInteger(genename, geneCount);
                    geneCount++;
                }
            }
        }
    }
    fclose(inFile);

    in = fopen(infile, "r");
    IntArray *SNPList;
    SNPList = new IntArray[geneCount];

    while (!feof(inFile))
    {
        buffer.ReadLine(inFile);
        tmp.Clear();
        tmp.AddTokens(buffer, "\t;");
        for (int i = 0; i < tmp.Length(); i++)
        {
            if (tmp[i].Find(function) != -1)
            {
                StringArray sub;
                sub.AddTokens(tmp[i], ":");
                String genename = sub[1];
                int gene = groupHash.Integer(genename);
                SNPList[gene].Push(tmp[0] + ":" + tmp[1]);
            }
        }
    }
    fclose(inFile);

    File *out = fopen(outfile, "wt");
    for (int i = 0; i < geneCount; i++)
    {
        fprintf(out, "%s\t", geneNames[i]);
        for (int m = 0; m < SNPList[i].Length(); m++)
        {
            fprintf(out, "%s\t", SNPList[i][m].c_str());
        }
        fprintf(out, "\n");
    }
    fclose(out);

    return 0;
}

