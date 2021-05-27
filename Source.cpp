#include "Header.h"

double sphere::volume() {
	double res = 4 * M_PI * pow(Rad, 3) / 3;
	return res;
}

double box::volume() {
	double res = pow(sqrt((pow(max[0] - min[0], 2) + pow(max[1] - min[1], 2) + pow(max[2] - min[2], 2)) / 3), 3);
	return res;
}

double tetra::volume() {
	double res = ((a[1][0] - a[0][0]) * (a[2][1] - a[0][1]) * (a[3][2] - a[0][2]) + (a[3][0] - a[0][0]) * (a[1][1] - a[0][1]) * (a[2][2] - a[0][2]) + (a[1][2] - a[0][2]) * (a[3][1] - a[0][1]) * (a[2][0] - a[0][0]) - (a[1][2] - a[0][2]) * (a[2][1] - a[0][1]) * (a[3][0] - a[0][0]) - (a[3][2] - a[0][2]) * (a[1][1] - a[0][1]) * (a[2][0] - a[0][0]) - (a[1][0] - a[0][0]) * (a[2][2] - a[0][2]) * (a[3][1] - a[0][1])) / 6;
	return res;
}