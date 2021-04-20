#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "Bin.h"

class GlobalPlacer
{
public:
    GlobalPlacer(Placement &placement);
    void place();
    void plotPlacementResult(const string outfilename, bool isPrompt = false);
    void randomPlace(); // An example of random placement implemented by TA
    void HPWLOptimizer();
    void DensityHandler();
    void BinHandler();
    void DensityCalculator();   

private:
    Placement &_placement;
    Bin *_bins;   
    void plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2);
};

#endif // GLOBALPLACER_H