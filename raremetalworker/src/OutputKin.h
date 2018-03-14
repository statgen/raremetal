#ifndef __OUTPUTKIN_H__
#define __OUTPUTKIN_H__

#include "Kinship.h"
#include "KinshipX.h"
#include "StringBasics.h"

class OutputKin
{
public:
    OutputKin()
    {};

    ~OutputKin()
    {};
    static bool outputKinX, outputKin, kinOnly;

    static void Write(Family &f, int before, int N);

    static void WriteX(Family &f, int before, int N);

    void Output(Pedigree &ped);
};

#endif
