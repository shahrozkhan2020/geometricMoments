#ifndef geometricMoments_H
#define geometricMoments_H
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <math.h> 
#include <cmath>
#include <fstream> 
#include <algorithm>

#include <utility> // pair
#include <stdexcept> // runtime_error
#include <sstream> // stringstream
#include <iomanip> 
#include "stl_reader.h"
#include <ppl.h>
#include <Windows.h>
#include <array>
#include <stdio.h>

using namespace std;
using namespace concurrency;

class geometricMoments
{
public:
	stl_reader::StlMesh <float, unsigned int> m_mesh;
public:
	geometricMoments(string stlGeometryPath);
	virtual ~geometricMoments();
public:
	double getMoment(double p, double q, double r, bool isCentered, bool isScaled);
	vector<pair<string, double>> getMomentVector(double order, bool isCentered, bool isScaled);
	
private:
	vector<vector<vector<double>>> momentMatrix(double p, double q, double r);
	double factorial(double n);
	void getVolumeProperties(vector<double>& centeriod, double& volume);
};

#endif 