#pragma once
#include <cmath>

extern const float M_PI = 3.14159265358979323846f;  // pi
extern const float MU_km = 3.986f * pow(10.0f,5.0f); // in km^3 / s^2
extern const float MU = 1;// in AU^3 / TU^2
#define sind(x) (sin((x) * M_PI / 180.0f))
#define cosd(x) (cos((x) * M_PI / 180.0f))

//#define asind(x) (asin((x) * M_PI / 180))
//#define acosd(x) (acos((x) * M_PI / 180))
