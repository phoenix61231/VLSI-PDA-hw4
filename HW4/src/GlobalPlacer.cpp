#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "Bin.h"
#include <math.h>

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement)
{
}

// Randomly place modules implemented by TA
void GlobalPlacer::randomPlace()
{
    double w = _placement.boundryRight() - _placement.boundryLeft();
    double h = _placement.boundryTop() - _placement.boundryBottom();
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        double wx = _placement.module(i).width(),
               hx = _placement.module(i).height();
        double px = (int)rand() % (int)(w - wx) + _placement.boundryLeft();
        double py = (int)rand() % (int)(h - hx) + _placement.boundryBottom();
        _placement.module(i).setPosition(px, py);
    }
}

void GlobalPlacer::HPWLOptimizer()
{
    HPWL_Function ef(_placement, 0); // require to define the object function and gradient function

    vector<double> x(2); // solution vector, size: num_blocks*2
                         // each 2 variables represent the X and Y dimensions of a block  
    int iter = 1;

    while(iter>0){
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        //_placement.module(i).setPosition(_placement.boundryLeft(), _placement.boundryBottom());

        //srand(time(NULL));
        //_placement.module(i).setPosition((int)rand() % (int)(W - _placement.module(i).width()) + _placement.boundryLeft(), (int)rand() % (int)(H - _placement.module(i).height()) + _placement.boundryBottom());
        
        ef.module_num = i;

        // initialize the solution vector
        x[0] = _placement.module(i).x();
        x[1] = _placement.module(i).y();

        NumericalOptimizer no(ef);
        no.setX(x);             // set initial solution
        no.setNumIteration(10); // user-specified parameter
        no.setStepSizeBound(3); // user-specified parameter
        no.solve();             // Conjugate Gradient solver

        double px = no.x(0);
        double py = no.x(1);

        if(px < _placement.boundryLeft()){
            px = _placement.boundryLeft();
        }
        else if(px > _placement.boundryRight()){
            px = _placement.boundryRight()-_placement.module(i).width();
        }

        if(py < _placement.boundryBottom()){
            py = _placement.boundryBottom();
        }
        else if(py > _placement.boundryTop()){
            py = _placement.boundryTop()-_placement.module(i).height();
        }

        cout << "Objective: " << no.objective() << endl;

        _placement.module(i).setPosition(px, py);
    }
    
    iter--;
    }
}

void GlobalPlacer::DensityHandler()
{
    double RATIO = 1.0;

    Bin *curr_bin = _bins, *tar_bin, *next_row, *srh_bin;
    double W = _placement.boundryRight() - _placement.boundryLeft(); 
    double H = _placement.boundryTop() - _placement.boundryBottom();

    double px, py;

    cout << "Density Handler" << endl;

    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        curr_bin = _bins;

        //////////  find the bin where the module locate    //////////
        while(_placement.module(i).x()>(curr_bin->_x+W/10) || _placement.module(i).x()<curr_bin->_x){
            curr_bin = curr_bin->_right;
        }

        while(_placement.module(i).y()>(curr_bin->_y+H/10) || _placement.module(i).y()<curr_bin->_y){
            curr_bin = curr_bin->_up;
        }
        //cout << curr_bin->_area << endl;

        if(curr_bin->_ratio>=RATIO){
            //////////  find a bin which density is small enough //////////
            /*tar_bin = _bins;
            next_row = _bins;
            while(((tar_bin->_area+_placement.module(i).area())/((W/10)*(H/10)))>RATIO){
                if(tar_bin->_right!=NULL){
                    tar_bin = tar_bin->_right;
                    continue;
                }

                if(next_row!=NULL){
                    tar_bin = next_row;
                    next_row = next_row->_up;
                }
                else break;
            }*/

            //////////  find a bin which density is small enough & nearest to current bin//////////
            srh_bin = _bins;
            next_row = _bins;
            tar_bin = _bins;

            while(1){
                if(((srh_bin->_area+_placement.module(i).area())/((W/10)*(H/10)))<RATIO){
                    double srh_d = pow(pow(curr_bin->_x-srh_bin->_x, 2)+pow(curr_bin->_y-srh_bin->_y, 2), 1/2);
                    double tar_d = pow(pow(curr_bin->_x-tar_bin->_x, 2)+pow(curr_bin->_y-tar_bin->_y, 2), 1/2);
                    
                    if(srh_d < tar_d){
                        tar_bin = srh_bin;
                    }
                }

                if(srh_bin->_right!=NULL){
                    srh_bin = srh_bin->_right;
                    continue;
                }

                if(next_row!=NULL){
                    srh_bin = next_row;
                    next_row = next_row->_up;
                }
                else break;
            }

            //////////  update density  //////////
            /*if((_placement.module(i).x()+_placement.module(i).width())<=(curr_bin->_x+W/10) && 
                (_placement.module(i).y()+_placement.module(i).height())<=(curr_bin->_y+H/10))
            {
                curr_bin->_area -= _placement.module(i).width()*_placement.module(i).height();
            }
            else if(((_placement.module(i).x()+_placement.module(i).width())<=(curr_bin->_x+W/10) && 
                (_placement.module(i).y()+_placement.module(i).height())>(curr_bin->_y+H/10)))
            {
                curr_bin->_area -= _placement.module(i).width()*((curr_bin->_y+H/10)-_placement.module(i).y());
                curr_bin->_up->_area -= _placement.module(i).width()*(_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y());
            }
            else if(((_placement.module(i).x()+_placement.module(i).width())>(curr_bin->_x+W/10) && 
                (_placement.module(i).y()+_placement.module(i).height())<=(curr_bin->_y+H/10)))
            {
                curr_bin->_area -= _placement.module(i).height()*((curr_bin->_x+W/10)-_placement.module(i).x());
                curr_bin->_right->_area -= _placement.module(i).height()*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
            }
            else{
                curr_bin->_area -= ((curr_bin->_x+W/10)-_placement.module(i).x())*((curr_bin->_y+H/10)-_placement.module(i).y());
                curr_bin->_up->_area -= ((curr_bin->_x+W/10)-_placement.module(i).x())*(_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y());
                curr_bin->_right->_area -= ((curr_bin->_y+H/10)-_placement.module(i).y())*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
                curr_bin->_right->_up->_area -= (_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y())*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
            }*/

            //////////  move (x_cor & y_cor) update density //////////
            px = (int)rand() % (int)(W/10 - _placement.module(i).width()) + tar_bin->_x;
            py = (int)rand() % (int)(H/10 - _placement.module(i).height()) + tar_bin->_y;

            //////////  update density  //////////
            curr_bin->_area -= _placement.module(i).width()*_placement.module(i).height();
            tar_bin->_area += _placement.module(i).width()*_placement.module(i).height();

            curr_bin->_ratio = curr_bin->_area/((W/10)*(H/10));
            tar_bin->_ratio = tar_bin->_area/((W/10)*(H/10));
            
            //px = _placement.module(i).x();
            //py = _placement.module(i).y();

            _placement.module(i).setPosition(px, py);            
        }     
    }
}

void GlobalPlacer::BinHandler()
{
    double W = _placement.boundryRight() - _placement.boundryLeft(); 
    double H = _placement.boundryTop() - _placement.boundryBottom();
    Bin *new_bin = NULL, *curr_row = NULL, *last_row = NULL, *last_bin = NULL;

    cout << "Bin Handler" << endl;

    for(double y_curr = _placement.boundryBottom(); y_curr<_placement.boundryTop(); y_curr+=H/10){
        for(double x_curr = _placement.boundryLeft(); x_curr<_placement.boundryRight(); x_curr+=W/10){
            new_bin = new Bin();
            if(x_curr == _placement.boundryLeft() && y_curr == _placement.boundryBottom()) _bins = new_bin;

            new_bin->_x = x_curr;
            new_bin->_y = y_curr;
            
            if(last_bin!=NULL){
                new_bin->_left = last_bin;
                last_bin->_right = new_bin;
            }

            if(last_row!=NULL){
                new_bin->_down = last_row;
                last_row->_up = new_bin;
                last_row = last_row->_right;                
            }

            if(curr_row == NULL){ curr_row = new_bin;}            

            last_bin = new_bin;

            //cout << "x : " << x_curr << " , y : " << y_curr << endl;
        }
        
        last_bin = NULL;
        last_row = curr_row;
        curr_row  = NULL;
    }
}

void GlobalPlacer::DensityCalculator()
{  
    double W = _placement.boundryRight() - _placement.boundryLeft(); 
    double H = _placement.boundryTop() - _placement.boundryBottom();
    double modules_area = 0, bins_area = 0;

    cout << "Density Calculator" << endl;

    Bin *curr_bin = _bins, *curr_row = _bins;
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        curr_bin = _bins;
        
        //////////  find the bin where the module locate    //////////
        while(_placement.module(i).x()>(curr_bin->_x+W/10) || _placement.module(i).x()<curr_bin->_x){
            curr_bin = curr_bin->_right;
        }

        while(_placement.module(i).y()>(curr_bin->_y+H/10) || _placement.module(i).y()<curr_bin->_y){
            curr_bin = curr_bin->_up;
        }    

        //////////  move modele, so whole module is in the bin  //////////
        double px = (int)rand() % (int)(W/10 - _placement.module(i).width()) + curr_bin->_x;
        double py = (int)rand() % (int)(H/10 - _placement.module(i).height()) + curr_bin->_y;

        curr_bin->_area += _placement.module(i).width()*_placement.module(i).height();
        _placement.module(i).setPosition(px, py);        

        //////////  calculate area distribution //////////
        /*if((_placement.module(i).x()+_placement.module(i).width())<=(curr_bin->_x+W/10) && 
            (_placement.module(i).y()+_placement.module(i).height())<=(curr_bin->_y+H/10))
        {
            curr_bin->_area += _placement.module(i).width()*_placement.module(i).height();
        }
        else if(((_placement.module(i).x()+_placement.module(i).width())<=(curr_bin->_x+W/10) && 
            (_placement.module(i).y()+_placement.module(i).height())>(curr_bin->_y+H/10)))
        {
            curr_bin->_area += _placement.module(i).width()*((curr_bin->_y+H/10)-_placement.module(i).y());
            curr_bin->_up->_area += _placement.module(i).width()*(_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y());
        }
        else if(((_placement.module(i).x()+_placement.module(i).width())>(curr_bin->_x+W/10) && 
            (_placement.module(i).y()+_placement.module(i).height())<=(curr_bin->_y+H/10)))
        {
            curr_bin->_area += _placement.module(i).height()*((curr_bin->_x+W/10)-_placement.module(i).x());
            curr_bin->_right->_area += _placement.module(i).height()*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
        }
        else{
            curr_bin->_area += ((curr_bin->_x+W/10)-_placement.module(i).x())*((curr_bin->_y+H/10)-_placement.module(i).y());
            curr_bin->_up->_area += ((curr_bin->_x+W/10)-_placement.module(i).x())*(_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y());
            curr_bin->_right->_area += ((curr_bin->_y+H/10)-_placement.module(i).y())*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
            curr_bin->_right->_up->_area += (_placement.module(i).height()-(curr_bin->_y+H/10)+_placement.module(i).y())*(_placement.module(i).width()-(curr_bin->_x+W/10)+_placement.module(i).x());
        }*/

        modules_area += _placement.module(i).height()*_placement.module(i).width();
    }


    curr_bin = _bins;
    //////////  Traversal   //////////
    cout << "Traverse Bins" << endl;
    while(curr_row!=NULL){
        while(curr_bin!=NULL){
            //cout << "x : " << curr_bin->_x << " ,y : " << curr_bin->_y << endl;
            bins_area += curr_bin->_area;
            curr_bin->_ratio = curr_bin->_area/((W/10)*(H/10));
            //cout << "ratio : " << curr_bin->_ratio << endl;

            curr_bin = curr_bin->_right;            
        }

        curr_row = curr_row->_up;
        curr_bin = curr_row;
    }   

    //cout << "module : " << modules_area << ", bin : " << bins_area << endl;
}

void GlobalPlacer::place()
{
    ///////////////////////////////////////////////////////////////////
    // The following example is only for analytical methods.
    // if you use other methods, you can skip and delete it directly.
    //////////////////////////////////////////////////////////////////   

    //////////  find initial placement  //////////
    //////////  use the conjugate optimizer //////////
    HPWLOptimizer();

    //////////  calculate bin density   //////////
    //////////  cut to bin  //////////
    BinHandler();
    DensityCalculator();

    //////////  don't use conjugate optimizer, use density  //////////
    DensityHandler();

    ////////////////////////////////////////////////////////////////

    // An example of random placement by TA. If you want to use it, please uncomment the folllwing 2 lines.
    // srand(time(NULL));
    // randomPlace();

    /* @@@ TODO 
	 * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
	 * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN 
	 * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
	 * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
	 * 7. Set the initial vector x in main(), set step size, set #iteration, and call the solver like above example
	 * */
}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt)
{
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt)
    {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd))
        {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2)
{
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}
