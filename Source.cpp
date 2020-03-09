#include <iostream>
#include <math.h>
#include "Header.h"
#include <cmath>
using namespace std;

// Determine planet centered position vector
void planetcentric(double (*vec)[3][1], double (*pcp)[3][1]) {
	for (int i = 0; i < 3; i++) {
		(*pcp)[i][0] = (*vec)[i][0];
	}
	//double Re = 6371.0088;//Re is the radius of earth, in this case we use the mean radius of earth in km
	double Re = 1.0; // Radius of earth in DU
	(*pcp)[2][0] = Re + (*pcp)[2][0];
	//cout << (*pcp)[2][0] << endl;
}


void transpose(double (*mat)[3][3], double (*d_transpose)[3][3]) {
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			(*d_transpose)[i][j] = (*mat)[j][i];
		}
	}
}


void matmult(double (*result)[3][1], double (*A)[3][3], double (*B)[3][1], int rowsA, int rowsB, int colsA, int colsB) {
	for (int i = 0; i < rowsA; i++) {
		for (int j = 0; j < colsB; j++) {
			(*result)[i][j] = 0;
			for (int k = 0; k < colsA; k++) {
				(*result)[i][j] += ((*A)[i][k]) * ((*B)[k][j]);
			}
		}
	}
}


void crossP(double (*cross)[3][1], double (*A)[3][1], double (*B)[3][1]) {
	// crossx = AyBz - AzBy
	(*cross)[0][0] = ( (*A)[1][0] * (*B)[2][0] ) - ( (*A)[2][0] * (*B)[1][0] );
	// crossy = AzBx - AxBz
	(*cross)[1][0] = ((*A)[2][0] * (*B)[0][0]) - ((*A)[0][0] * (*B)[2][0]);
	// crossz = AxBy - AyBx
	(*cross)[2][0] = ((*A)[0][0] * (*B)[1][0]) - ((*A)[1][0] * (*B)[0][0]);
}

double magnitude (double (*vec)[3][1]) {
	double x1 = 0;
	x1 = pow(((*vec)[0][0]), 2);
	double x2 = 0;
	x2 = pow(((*vec)[1][0]), 2);
	double x3 = 0;
	x3 = pow(((*vec)[2][0]), 2);

	double mag = 0;
	mag = sqrt( x1 + x2 + x3 );
	//cout << endl << mag << endl;
	return mag;
}

void scalarMult (double (*newVec)[3][1], double(*vec)[3][1], double scalar) {
	for (int i = 0; i < 3; i++) {
		(*newVec)[i][0] += (*vec)[i][0] * scalar;
	}
}

double dotP (double(*vec1)[3][1], double(*vec2)[3][1] ) {
	double scalar = 0;
	for (int i = 0; i < 3; i++) {
		scalar += (*vec1)[i][0] * (*vec2)[i][0];
	}
	return scalar;
}


int main() {
	// First we define the range and range rate measurements in the topocentric horizon coordinate system, lat and lon on radar station
	double rangeXh, rangeYh, rangeZh;
	double rateXh, rateYh, rateZh;

	double lat; // latitude of radar site
	double lambdaE; // This is the geographic longitude of the radar site
	double thetaGnot; // Greenwich mean sidereal time at known time, t0 (based off J-2000 epoch)
	double days; // This is the number of days after the known time t0

	// Take input for range, range rate, lat and lon of radar site, initial greenwich sidereal time and days past greenwich sidereal time

	/*cout << "Enter range in the topocentric horizon frame \nXh: ";
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
	cin >> days;*/

	rangeXh = -0.033280;
	rangeYh = 0.041224;
	rangeZh = 0.048642;
	rateXh = 0.663693;
	rateYh = 0.663850;
	rateZh = -0.003025;
	lat = 28.48417;
	lambdaE = -80.5725;
	thetaGnot = 100.36245;
	days = 124.3333;

	double rangeRate[3][1] = { { rateXh }, {rateYh}, {rateZh} };
	double range[3][1] = { {rangeXh}, {rangeYh}, {rangeZh} };

	int sizeVec = sizeof(range) / sizeof(range[0][0]); // Size of array is given by the (number of bytes) / (bytes in one index)

	/*cout << "range vector is: ";
	for (int i = 0; i < sizeVec; i++)
	{
		cout << range[i][0] << " ";
	}
	cout << "\n";
	cout << "range rate vector is: ";
	for (int i = 0; i < sizeVec; i++)
	{
		cout << rangeRate[i][0] << ' ';
	}
	cout << "\n\n";*/

	double pcp[3][1] = { {0}, {0}, {0} };
	planetcentric(&range, &pcp);
	//cout << *pcp[0] << "\n" << *pcp[1] << "\n" << *pcp[2];
	/*double newRange[3][1] = { {ptr[0]}, {ptr[1]}, {ptr[2]} };
	cout << newRange[0][0] << " " << newRange[1][0] << " " << newRange[2][0] << endl;*/

	double theta;
	double thetaG = thetaGnot + (1.00273779093 * 360 * days);
	int thetaGint = trunc(thetaG);
	theta = thetaGint % 360;
	double diff = thetaG - thetaGint;
	theta = theta + diff;
	theta = theta + lambdaE;
	//cout << theta << endl;

	//double d[3][3]; // this is the D matrix that transforms Earth centered inertial frame to topocentric horizon frame so we need the inverse since we know the topocentric horizon frame
	double r11 = sind(lat) * cosd(theta);
	double r12 = sind(lat) * sind(theta);
	double r13 = -1 * cosd(lat);
	double r21 = -1 * sind(theta);
	double r22 = cosd(theta);
	double r23 = 0;
	double r31 = cosd(lat) * cosd(theta);
	double r32 = cosd(lat) * sind(theta);
	double r33 = sind(lat);
	//cout << r11 << " " << r12 << " " << r13 << endl << r21 << " " << r22 << " " << r23 << endl << r31 << " " << r32 << " " << r33 << "\n\n";
	double d[3][3] = { { r11, r12, r13 },{ r21, r22, r23 }, { r31, r32, r33 } };


	// We need the inverse of D to get X,Y,Z cooridinate frame and D inverse = D transpose in this case
	double d_transpose[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
	transpose(&d, &d_transpose);

	/*printf("Result matrix is \n");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			cout << d_transpose[i][j] << " ";
		}
		cout << endl;
	}*/

	double range_pc[3][1] = { {0},{0},{0} };
	int rowsA = 3, rowsB = 3, colsA = 3, colsB = 1;
	matmult(&range_pc, &d_transpose, &pcp, rowsA, rowsB, colsA, colsB);

	/*cout << endl;
	for (int i = 0; i < 3; i++)
	{
		int j = 0;
		cout << range_pc[i][j] << " ";
		cout << endl;
	}*/

	double rate_pc[3][1] = { {0}, {0}, {0} };
	matmult(&rate_pc, &d_transpose, &rangeRate, rowsA, rowsB, colsA, colsB);

	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		int j = 0;
		cout << rate_pc[i][j];
		cout << endl;
	}*/

	// Calculate velocity in the planet centered frame by including coriolis acceleration
	// Vel = rangeRate + (omega X range), must make a cross product function
	double w[3][1] = { {0}, {0}, {0.05883} }; // Units of rad / TU
	double cross[3][1] = { {0}, {0}, {0} };
	crossP(&cross, &w, &range_pc);

	//double vel[3][1] = { {*cross[0][0]}, {*cross[1][0]}, {*cross[2][0]} };
	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		int j = 0;
		cout << cross[i][j];
		cout << endl;
	}*/

	double vel[3][1] = { {0}, {0}, {0} };
	for (int i = 0; i < 3; i++) {
		int j = 0;
		(vel)[i][j] = (rate_pc)[i][j] + (cross)[i][j];
	}
	
	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		cout << vel[i][0] << endl;
	}*/

	// Now we must compute the three fundamental vectors
	// angular momentum = r X V
	double h[3][1] = { {0}, {0}, {0} };
	double h_mag = 0;

	crossP(&h, &range_pc, &vel);

	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		cout << h[i][0] << endl;
	}*/
	h_mag = magnitude(&h);
	//cout << endl << h_mag;


	// calculate n = Z X h and magnitude of n
	double z[3][1] = { {0}, {0}, {1} };
	double n[3][1] = { {0}, {0}, {0} };
	double n_mag = 0;

	crossP(&n, &z, &h);

	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		cout << n[i][0] << endl;
	}*/
	n_mag = magnitude(&n);
	//cout << endl << n_mag << endl;

	// calculate e vector and magnitude
	double e[3][1] = { {0}, {0}, {0} };
	double e_mag = 0;

	double rdotv = dotP(&range_pc, &vel);
	//cout << endl << rdotv << endl;
	double v_mag = magnitude(&vel);

	double r = magnitude(&range_pc);
	for (int i = 0; i < 3; i++) {
		e[i][0] = (1 / MU) * (((pow(v_mag, 2) - (MU / r)) * (range_pc)[i][0]) - (rdotv * (vel)[i][0]));
		//cout << endl << e[i][0] << endl;
	}
	e_mag = magnitude(&e);

	//cout << endl << e_mag << endl;

	// Now we calculate the orbital elements

	double p = pow(h_mag, 2) / MU; //cout << endl << p;

	double a = p / (1 - pow(e_mag, 2));

	double ndote = dotP(&n, &e);
	double loweromega = acos(ndote / (n_mag * e_mag));
	loweromega = loweromega * (180 / M_PI);
	
	double upperomega = acos(((n[0][0]) / n_mag));
	upperomega = upperomega * (180 / M_PI);

	double i = acos((h[2][0]) / h_mag);
	i = i * (180 / M_PI);

	double edotr = dotP(&e, &range_pc);
	double nu = acos( edotr / (e_mag * r));
	nu = nu * (180 / M_PI);

	cout << "The position vector in ECI frame: " << range_pc[0][0] << " X + " << range_pc[1][0] << " Y + " << range_pc[2][0] << " Z  DU" << endl;
	cout << "THe velocity vector in ECI frome: " << vel[0][0] << " X + " << vel[1][0] << " Y + " << vel[2][0] << " Z  DU/TU" << "\n\n";

	cout << "The Semi-major axis, a: " << a << " DU" << endl;
	cout << "The Eccentricity, e: " << e_mag << endl;
	cout << "The Inclination i: " << i << endl;
	cout << "The Argument of Perigee, lowercase omega: " << loweromega << endl;
	cout << "The Right ascention of the Ascending Node, capital omega: " << upperomega << endl;
	cout << "The True Anomaly, Nu: " << nu << endl;
	//cout << endl << "The orbital element p: " << p << endl << "The orbital element a: " << a << endl << "The orbital element \u03C9: " << loweromega << endl << "The orbital element \u03A9:" << upperomega << endl << "The orbital element i: " << i << endl << "The orbital element \u03BD:" << nu << endl;
}