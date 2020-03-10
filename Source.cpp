#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;

// Determine planet centered position vector
void planetcentric(float (*vec)[3][1], float (*pcp)[3][1]) {
	for (int i = 0; i < 3; i++) {
		(*pcp)[i][0] = (*vec)[i][0];
	}
	//float Re = 6371.0088;//Re is the radius of earth, in this case we use the mean radius of earth in km
	float Re = 1.0; // Radius of earth in DU
	(*pcp)[2][0] = Re + (*pcp)[2][0];
}


void transpose(float (*mat)[3][3], float (*d_transpose)[3][3]) {
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			(*d_transpose)[i][j] = (*mat)[j][i];
		}
	}
}


void matmult(float (*result)[3][1], float (*A)[3][3], float (*B)[3][1], int rowsA, int rowsB, int colsA, int colsB) {
	for (int i = 0; i < rowsA; i++) {
		for (int j = 0; j < colsB; j++) {
			(*result)[i][j] = 0;
			for (int k = 0; k < colsA; k++) {
				(*result)[i][j] += ((*A)[i][k]) * ((*B)[k][j]);
			}
		}
	}
}

void crossP(float (*cross)[3][1], float (*A)[3][1], float (*B)[3][1]) {
	// crossx = AyBz - AzBy
	(*cross)[0][0] = ( (*A)[1][0] * (*B)[2][0] ) - ( (*A)[2][0] * (*B)[1][0] );
	// crossy = AzBx - AxBz
	(*cross)[1][0] = ((*A)[2][0] * (*B)[0][0]) - ((*A)[0][0] * (*B)[2][0]);
	// crossz = AxBy - AyBx
	(*cross)[2][0] = ((*A)[0][0] * (*B)[1][0]) - ((*A)[1][0] * (*B)[0][0]);
}

float magnitude (float (*vec)[3][1]) {
	float x1 = 0;
	x1 = pow(((*vec)[0][0]), 2);
	float x2 = 0;
	x2 = pow(((*vec)[1][0]), 2);
	float x3 = 0;
	x3 = pow(((*vec)[2][0]), 2);

	float mag = 0;
	mag = sqrt( x1 + x2 + x3 );
	//cout << endl << mag << endl;
	return mag;
}

void scalarMult (float (*newVec)[3][1], float(*vec)[3][1], float scalar) {
	for (int i = 0; i < 3; i++) {
		(*newVec)[i][0] += (*vec)[i][0] * scalar;
	}
}

float dotP (float(*vec1)[3][1], float(*vec2)[3][1] ) {
	float scalar = 0;
	for (int i = 0; i < 3; i++) {
		scalar += (*vec1)[i][0] * (*vec2)[i][0];
	}
	return scalar;
}


int main() {
	// First we define the range and range rate measurements in the topocentric horizon coordinate system, lat and lon on radar station
	float rangeXh, rangeYh, rangeZh;
	float rateXh, rateYh, rateZh;

	float lat; // latitude of radar site
	float lambdaE; // This is the geographic longitude of the radar site
	float thetaGnot; // Greenwich mean sidereal time at known time, t0 (based off J-2000 epoch)
	float days; // This is the number of days after the known time t0

	// Take input for range, range rate, lat and lon of radar site, initial greenwich sidereal time and days past greenwich sidereal time
	cout << "Enter range in the topocentric horizon frame in DU \nXh: ";
	cin >> rangeXh;
	cout << "Yh: ";
	cin >> rangeYh;
	cout << "Zh: ";
	cin >> rangeZh;
	//float range[3] = { rangeXh, rangeYh, rangeZh };

	cout << "Enter range rate in the topocentric horizon frame in DU/TU \nXh: ";
	cin >> rateXh;
	cout << "Yh: ";
	cin >> rateYh;
	cout << "Zh: ";
	cin >> rateZh;
	//float rangeRate[3] = { rateXh, rateYh, rateZh };

	cout << "Enter latitude of the radar site in degrees (+'ve North): ";
	cin >> lat;
	cout << "Enter longitude of the radar site in degrees (+'ve East): ";
	cin >> lambdaE;
	cout << "Enter the original Greenwich sidereal time in degrees: ";
	cin >> thetaGnot;
	cout << "Enter the number of days after past this original Greenwich sidereal time in decimal form: ";
	cin >> days;

	/*rangeXh = -0.033280f;
	rangeYh = 0.041224f;
	rangeZh = 0.048642f;
	rateXh = 0.663693f;
	rateYh = 0.663850f;
	rateZh = -0.003025f;
	lat = 28.48417f;
	lambdaE = -80.5725f;
	thetaGnot = 100.36245f;
	days = 124.3333f;*/

	float rangeRate[3][1] = { { rateXh }, {rateYh}, {rateZh} };
	float range[3][1] = { {rangeXh}, {rangeYh}, {rangeZh} };

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

	float pcp[3][1] = { {0}, {0}, {0} };
	planetcentric(&range, &pcp);

	float theta;
	float thetaG = thetaGnot + (1.00273779093f * 360.0f * days);
	int thetaGint = int(trunc(thetaG));
	theta = float(int(thetaGint % 360));
	float diff = thetaG - thetaGint;
	theta = theta + diff;
	theta = theta + lambdaE;
	//cout << theta << endl;

	//float d[3][3]; // this is the D matrix that transforms Earth centered inertial frame to topocentric horizon frame so we need the inverse since we know the topocentric horizon frame
	float r11 = float(sind(lat) * cosd(theta));
	float r12 = float(sind(lat) * sind(theta));
	float r13 = float(-1 * cosd(lat));
	float r21 = float(-1 * sind(theta));
	float r22 = float(cosd(theta));
	float r23 = 0.0f;
	float r31 = float(cosd(lat) * cosd(theta));
	float r32 = float(cosd(lat) * sind(theta));
	float r33 = float(sind(lat));
	//cout << r11 << " " << r12 << " " << r13 << endl << r21 << " " << r22 << " " << r23 << endl << r31 << " " << r32 << " " << r33 << "\n\n";
	float d[3][3] = { { r11, r12, r13 },{ r21, r22, r23 }, { r31, r32, r33 } };


	// We need the inverse of D to get X,Y,Z cooridinate frame and D inverse = D transpose in this case
	float d_transpose[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
	transpose(&d, &d_transpose);

	/*printf("Result matrix is \n");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			cout << d_transpose[i][j] << " ";
		}
		cout << endl;
	}*/

	float range_pc[3][1] = { {0},{0},{0} };
	int rowsA = 3, rowsB = 3, colsA = 3, colsB = 1;
	matmult(&range_pc, &d_transpose, &pcp, rowsA, rowsB, colsA, colsB);

	/*cout << endl;
	for (int i = 0; i < 3; i++)
	{
		int j = 0;
		cout << range_pc[i][j] << " ";
		cout << endl;
	}*/

	float rate_pc[3][1] = { {0}, {0}, {0} };
	matmult(&rate_pc, &d_transpose, &rangeRate, rowsA, rowsB, colsA, colsB);

	// Calculate velocity in the planet centered frame by including coriolis acceleration
	// Vel = rangeRate + (omega X range), must make a cross product function
	float w[3][1] = { {0}, {0}, {0.05883f} }; // Units of rad / TU
	float cross[3][1] = { {0}, {0}, {0} };
	crossP(&cross, &w, &range_pc);

	float vel[3][1] = { {0}, {0}, {0} };
	for (int i = 0; i < 3; i++) {
		int j = 0;
		(vel)[i][j] = (rate_pc)[i][j] + (cross)[i][j];
	}

	// Now we must compute the three fundamental vectors
	// angular momentum = r X V
	float h[3][1] = { {0}, {0}, {0} };
	float h_mag = 0;

	crossP(&h, &range_pc, &vel);

	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		cout << h[i][0] << endl;
	}*/
	h_mag = magnitude(&h);
	//cout << endl << h_mag;


	// calculate n = Z X h and magnitude of n
	float z[3][1] = { {0}, {0}, {1} };
	float n[3][1] = { {0}, {0}, {0} };
	float n_mag = 0;

	crossP(&n, &z, &h);

	/*cout << endl;
	for (int i = 0; i < 3; i++) {
		cout << n[i][0] << endl;
	}*/
	n_mag = magnitude(&n);
	//cout << endl << n_mag << endl;

	// calculate e vector and magnitude
	float e[3][1] = { {0}, {0}, {0} };
	float e_mag = 0;

	float rdotv = dotP(&range_pc, &vel);
	//cout << endl << rdotv << endl;
	float v_mag = magnitude(&vel);

	float r = magnitude(&range_pc);
	for (int i = 0; i < 3; i++) {
		e[i][0] = (1 / MU) * (((pow(v_mag, 2) - (MU / r)) * (range_pc)[i][0]) - (rdotv * (vel)[i][0]));
		//cout << endl << e[i][0] << endl;
	}
	e_mag = magnitude(&e);

	//cout << endl << e_mag << endl;

	// Now we calculate the orbital elements

	float p = pow(h_mag, 2) / MU; //cout << endl << p;

	float a = p / (1 - pow(e_mag, 2));

	float ndote = dotP(&n, &e);
	float loweromega = acos(ndote / (n_mag * e_mag));
	loweromega = loweromega * (180 / M_PI);
	
	float upperomega = acos(((n[0][0]) / n_mag));
	upperomega = upperomega * (180 / M_PI);

	float i = acos((h[2][0]) / h_mag);
	i = i * (180 / M_PI);

	float edotr = dotP(&e, &range_pc);
	float nu = acos( edotr / (e_mag * r));
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

	// Convert semi-major axis to Km: 1 DU = 6378.1 km
	float a_km = a * 6378.1f;

	// Calculate period in seconds
	//float period = (2 * M_PI) * sqrt((pow(a_km, 3)) / (MU_km));

	float t_pass = 50.5 * 60 * 60;

	float passages = 0;


	//float e = 0.8522;
	float Ei = 0.5; //This is a trial value for eccentric anomaly. Rule of thumb: Start with (Ei = pi) if e is large ( > 0.5). ** In Radians**
	//float mu = 42820; //Units of Km^3 / s^2
	//float a = 25823.16; //Units of Km
	float period = (2 * M_PI) * sqrt((pow(a_km,3))/(MU_km)); //Units of s. 35 Hours

	//float ToF = period - (float(10) * float(60)); // ToF (T - 20 mins from periapsis = Period - 10 mins). Units of s

	//float ToF = 50.5 * 3600;

	//nu = 66;

	//float ToF = 1565.88f; //This is the ToF for the example ploblem provided
	float ToF = 30.0f * 60.0f; //ToF of 30 minutes past initial measurement

	float n_meanmot = sqrt( (MU_km / (pow(a_km, 3))) );
	//float M = ToF * n_meanmot; //Calculated as ToF * n where n = sqrt(mu / a^3). n has units of Radians / s

	float E0 = acos( (e_mag + cosd(nu)) / (1 + (e_mag * (cosd(nu))) ));
	float M0 = E0 - (e_mag * sin(E0));
	float M = (ToF * n_meanmot) + M0;

	float Mi = 0; /////////////////////

	float M_Mi_Diff = 1;
	float Eip1 = 0;
	float epsilon = 1.0f * pow(10.0f, -7.0f);
	float dMdE = 0;
	float iteration = 0; //Keeps track of iteration number
	float Ei_temp = 0;
	cout << "\n";

	while (abs(M_Mi_Diff) > epsilon) {
		Ei_temp = Ei;
		// Step 1: Compute Mi
		Mi = Ei - e_mag * sin(Ei);
		M_Mi_Diff = M - Mi;
		// Step 2: Compute Eip1. i.e. Ei+1
		dMdE = 1 - (e_mag * cos(Ei));
		Ei = Ei + (M_Mi_Diff / (dMdE)); // 1 - e * cos Ei is dM/dE evaluated at Ei

		// Step 3: End when M-Mi is small. Included floato while loop condition. Answer is Ei+1
		cout << "Ei = " << Ei_temp << " Radians in iteration " << iteration << ", dMdE = " << dMdE << ", M diff = " << M_Mi_Diff << ", Ei+1 = " << Ei << "\n";
		iteration++;
	}

	float nu_iter = acos((cos(Ei) - e_mag) / (1 - (e_mag * cos(Ei))));

	if (Ei > M_PI) {
		nu_iter = (2 * M_PI) - nu_iter;
	}
	nu_iter = nu_iter * (180.0f / M_PI); //Convert to degrees
	cout << "\nNew Nu is: " << nu_iter << " degrees" << "\n";
	cout << "\nCompleted in " << iteration << " iterations\n";


	// Calculate position and velocity in perifocal-eccentricity frame (we will use equation 49)
	float rw[3][1] = { {float(r * cosd(nu_iter))}, {float(r * sind(nu_iter))}, {0} };

	float vw[3][1] = { {float((sqrt(MU / p)) * (-1) * sind(nu_iter))}, {float((sqrt(MU / p)) * (e_mag + cosd(nu_iter)))}, {0.0f} };

	cout << "The position vector in Perifocal-eccentricity frame: " << rw[0][0] << " Xw + " << rw[1][0] << " Yw DU" << endl;
	cout << "THe velocity vector in Perifocal-eccentricity frome: " << vw[0][0] << " Xw + " << vw[1][0] << " Yw DU/TU" << "\n\n";

	// Finally we need to convert this into planet-centered coordinate frame with transformation matrix R

	float R[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

	R[0][0] = float((cosd(upperomega) * cosd(loweromega)) - (sind(upperomega) * sind(loweromega) * cosd(i)));
	R[0][1] = float((-1 * cosd(upperomega) * sind(loweromega)) - (sind(upperomega) * cosd(loweromega) * cosd(i)));
	R[0][2] = float(sind(upperomega) * sind(i));
	R[1][0] = float((sind(upperomega) * cosd(loweromega)) + (cosd(upperomega) * sind(loweromega) * cosd(i)));
	R[1][1] = float((-1 * sind(upperomega) * sind(loweromega)) + (cosd(upperomega) * cosd(loweromega) * cosd(i)));
	R[1][2] = float(-1 * cosd(upperomega) * sind(i));
	R[2][0] = float(sind(loweromega) * sind(i));
	R[2][1] = float(cosd(loweromega) * sind(i));
	R[2][2] = float(cosd(i));

	float r_new[3][1] = { {0}, {0}, {0} };
	float v_new[3][1] = { {0}, {0}, {0} };

	matmult(&r_new, &R, &rw, 3, 3, 3, 1);
	matmult(&v_new, &R, &vw, 3, 3, 3, 1);

	cout << "The new position vector in ECI frame: " << r_new[0][0] << " X + " << r_new[1][0] << " Y + " << r_new[2][0] << " Z  DU" << endl;
	cout << "THe new velocity vector in ECI frome: " << v_new[0][0] << " X + " << v_new[1][0] << " Y + " << v_new[2][0] << " Z  DU/TU" << "\n\n";
}