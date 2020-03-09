#pragma once
#include <cmath>

extern const double M_PI = 3.14159265358979323846;  // pi
//extern const double MU = 3.986 * pow(10,5); // in km^3 / s^2
extern const double MU = 1;// in AU^3 / TU^2
#define sind(x) (sin((x) * M_PI / 180))
#define cosd(x) (cos((x) * M_PI / 180))

//#define asind(x) (asin((x) * M_PI / 180))
//#define acosd(x) (acos((x) * M_PI / 180))
