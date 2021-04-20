#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction();

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
};


class HPWL_Function : public NumericalOptimizerInterface
{
public:
    HPWL_Function(Placement &placement, int module_num);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

    Placement &_placement;
    int module_num;
};

#endif // EXAMPLEFUNCTION_H
