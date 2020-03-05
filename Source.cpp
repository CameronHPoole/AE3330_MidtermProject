#include <iostream>
using namespace std;

// Determine planet centered position vector
double* planetcentric(double* vec) {
	double pcp[3] = { vec[0], vec[1], vec[2] };
	double Re = 6371.0088;//Re is the radius of earth, in this case we use the mean radius of earth in km
	//cout << pcp[0];
	pcp[2] = Re + pcp[2];
	return pcp;
}

int main() {
	// First we define the range and range rate measurements in the topocentric horizon coordinate system, lat and lon on radar station
	double rangeXh;
	double rangeYh;
	double rangeZh;

	double rateXh;
	double rateYh;
	double rateZh;

	double lat; // latitude of radar site
	double lambdaE; // This is the geographic longitude of the radar site
	double thetaGnot; // Greenwich mean sidereal time at known time, t0 (based off J-2000 epoch)
	double days; // This is the number of days after the known time t0

	// Take input for range, range rate, lat and lon of radar site, initial greenwich sidereal time and days past greenwich sidereal time
	cout << "Enter range in the topocentric horizon frame \nXh: ";
	cin >> rangeXh;
	cout << "Yh: ";
	cin >> rangeYh;
	cout << "Zh: ";
	cin >> rangeZh;
	double range[3] = { rangeXh, rangeYh, rangeZh };

	cout << "Enter range rate in the topocentric horizon frame \nXh: ";
	cin >> rateXh;
	cout << "Yh: ";
	cin >> rateYh;
	cout << "Zh: ";
	cin >> rateZh;
	double rangeRate[3] = { rateXh, rateYh, rateZh };

	cout << "Enter latitude of the radar site in degrees: ";
	cin >> lat;
	cout << "Enter longitude of the radar site in degrees: ";
	cin >> lambdaE;
	cout << "Enter the original Greenwich sidereal time: ";
	cin >> thetaGnot;
	cout << "Enter the number of days after past this original Greenwich sidereal time: ";
	cin >> days;

	int sizeVec = sizeof(range) / sizeof(range[0]); // Size of array is given by the (number of bytes) / (bytes in one index)

	cout << "range vector is: ";
	for (int i = 0; i < sizeVec; i++)
	{
		cout << range[i] << " ";
	}
	cout << "\n";
	cout << "range rate vector is: ";
	for (int i = 0; i < sizeVec; i++)
	{
		cout << rangeRate[i] << ' ';
	}
	cout << "\n";

	double* ptr = planetcentric(range);
	//cout << ptr[0] << "\n" << ptr[1] << "\n" << ptr[2];
	
	double theta = thetaGnot + (1.002738 * 360 * days) + lambdaE;
	double d[3][3]; // this is the D matrix that transforms Earth centered inertial frame to topocentric horizon frame so we need the inverse since we know the topocentric horizon frame

}