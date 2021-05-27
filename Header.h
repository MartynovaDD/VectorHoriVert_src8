#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include "CImg.h"
using namespace cimg_library;
using namespace std;

class figure {
protected:
	double RGB[3];
	double XYZ[3];
public:
	virtual double volume() = 0;
	figure() {
		RGB[0] = 0; RGB[1] = 0; RGB[2] = 0;
		XYZ[0] = 0; XYZ[1] = 0; XYZ[2] = 0;
	}
	void SetColour(double clr[3]) { 
		for (int i = 0; i < 3; ++i) { 
		RGB[i] = clr[i]; 
		} 
	}
	friend class CRead;
};

class sphere : public figure {
private:
	double Rad;
public:
	sphere() { Rad = 0; }
	sphere(double x, double y, double z, double Rd) { 
		XYZ[0] = x; XYZ[1] = y; XYZ[2] = z; 
		Rad = Rd; 
	}
	double volume() override;
	friend class CRead;
};

class box : public figure {
private:
	double max[3];
	double min[3];
public:
	box() { 
		for (int i = 0; i < 3; ++i) { 
			min[i] = 0; max[i] = 0; 
		} 
	}
	box(double Min[3], double Max[3]) {
		for (int i = 0; i < 3; ++i) {
			min[i] = Min[i];
			max[i] = Max[i];
			XYZ[i] = (min[i] + max[i]) / 2;
		}
	}
	double volume() override;
	friend class CRead;
};

class tetra : public figure {
private:
	double a[4][3];
public:
	tetra() {
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				a[i][j] = 0;
			}
		}
	}
	tetra(double tetr[4][3]) {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 3; ++j) {
				a[i][j] = tetr[i][j];
			}
		}
		for (int i = 0; i < 3; ++i) {
			XYZ[i] = (tetr[0][i] + tetr[1][i] + tetr[2][i] + tetr[3][i]) / 4;
		}
	}
	double volume() override;
	friend class CRead;
};