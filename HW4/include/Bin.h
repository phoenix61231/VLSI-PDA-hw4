#ifndef BIN_H
#define BIN_H

using namespace std;

class Bin
{
public:
    Bin(){
        _x = 0;
        _y = 0;
        _ratio = 0;
        _area = 0;
        _up = NULL;
        _down = NULL;
        _left = NULL;
        _right = NULL;
    }

    double _x, _y;
    double _ratio;
    double _area;
    
    Bin *_up, *_down, *_left, *_right;    
};

#endif // BIN_H