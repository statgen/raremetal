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

    /**
     * `outputKin` controls whether to save the estimated GRM that is calculated (if --kinGeno is specified)
     * Saves the estimated GMR to a file named: yourprefix.Empirical.Kinship.gz
     * If --vcX is also issued in the command line, then a separate file named yourprefix.Empirical.KinshipX.gz will be
     *  generated where the GRM of chromosome X is saved.
     *
     * `kinOnly` allows users to estimate kinship matrix without any association analysis of any traits included in
     *  the data set. If specified, it will automatically cause the output to be saved as above.
     * To also estimate chromosome X kinship, --vcX option should be added in command line.
     */
    static bool outputKinX, outputKin, kinOnly;  // TODO: outputKinX may be totally unused/ ignored

    static void Write(Family &f, int before, int N);

    static void WriteX(Family &f, int before, int N);

    void Output(Pedigree &ped);
};

#endif
