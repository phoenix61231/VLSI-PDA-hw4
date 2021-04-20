#include "ExampleFunction.h"
#include <math.h>
#include "Pin.h"

ExampleFunction::ExampleFunction()
{
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function
    g[0] = 6 * x[0] + 2 * x[1];                                  // gradient function of X
    g[1] = 2 * x[0] + 4 * x[1];                                  // gradient function of Y
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function
}

unsigned ExampleFunction::dimension()
{
    return 2; // num_blocks*2
    // each two dimension represent the X and Y dimensions of each block
}

//////////  HPWL Function   //////////
HPWL_Function::HPWL_Function(Placement &placement, int module_num)
    : _placement(placement), module_num(module_num)
{
}

void HPWL_Function::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    f = 0;
    g[0] = 0;
    g[1] = 0;

    for(size_t i=0; i<_placement.module(module_num).numPins(); ++i){
        f += pow(_placement.module(module_num).pin(i).x()-x[0], 2);
        f += pow(_placement.module(module_num).pin(i).y()-x[1], 2);
    }  

    for(size_t i=0; i<_placement.module(module_num).numPins(); ++i){
        g[0] += 2*(_placement.module(module_num).pin(i).x()-x[0]);
    }  

    for(size_t i=0; i<_placement.module(module_num).numPins(); ++i){
        g[1] += 2*(_placement.module(module_num).pin(i).y()-x[1]);
    } 
}

void HPWL_Function::evaluateF(const vector<double> &x, double &f)
{
    f = 0;
    
    for(size_t i=0; i<_placement.module(module_num).numPins(); ++i){
        f += pow(_placement.module(module_num).pin(i).x()-x[0], 2);
        f += pow(_placement.module(module_num).pin(i).y()-x[1], 2);
    }     
}

unsigned HPWL_Function::dimension()
{
    return 2; // num_blocks*2
    // each two dimension represent the X and Y dimensions of each block
}