#include "geometricMoments.h"


geometricMoments::geometricMoments(string stlGeometryPath) {
	stl_reader::StlMesh <float, unsigned int> mesh(stlGeometryPath);
	m_mesh = mesh;
}

geometricMoments::~geometricMoments() {
}

double geometricMoments::getMoment(double p, double q, double r, bool isCentered, bool isScaled) {
	double s = p + q + r;
	double moment = 0;
	vector<vector<vector<double>>> s_out = momentMatrix(p, q, r);
	
	double m_000 = 0, m_010 = 0, m_100 = 0, m_001 = 0;
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, det;
	vector<double> centeriod; double volume;

	// If we want only scaled volume
	if (isScaled && p == 0 && q == 0 && r == 0) {
		return 1.0;
	}
	// we want centered first order
	if (isCentered && ((p == 1 && q == 0 && r == 0) || (p == 0 && q == 1 && r == 0) || (p == 0 && q == 0 && r == 1))) {
		return 0.0;
	}

	getVolumeProperties(centeriod, volume);


	// if we ask only volune 
	if (p == 0 && q == 0 && r == 0) {
		return volume;
	}
	// if we ask only first order because we are already calculating it - we do it only for nonscale case
	if (!isScaled) {
		if (p == 1 && q == 0 && r == 0) {
			return centeriod[0] * volume;
		}
		if (p == 0 && q == 1 && r == 0) {
			return centeriod[1] * volume;
		}
		if (p == 0 && q == 0 && r == 1) {
			return centeriod[2] * volume;
		}
	}
	if (isScaled) {
		if (p == 1 && q == 0 && r == 0) {
			return ((centeriod[0] * volume) / double((pow(volume, (1 + ((s) / 3.0))))));
		}
		if (p == 0 && q == 1 && r == 0) {
			return ((centeriod[1] * volume) / double((pow(volume, (1 + ((s) / 3.0))))));
		}
		if (p == 0 && q == 0 && r == 1) {
			return ((centeriod[2] * volume) / double((pow(volume, (1 + ((s) / 3.0))))));
		}
	}
	
	for (size_t itri = 0; itri < m_mesh.num_tris(); ++itri) {
		vector<vector<double>> A;
		const float* c;
		vector<double> A_temp;

		c = m_mesh.tri_corner_coords(itri, 0);
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		c = m_mesh.tri_corner_coords(itri, 1);
		A_temp.clear();
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		c = m_mesh.tri_corner_coords(itri, 2);
		A_temp.clear();
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		if (isCentered) {
			a1 = A[0][0] - centeriod[0]; a2 = A[0][1] - centeriod[1]; a3 = A[0][2] - centeriod[2];
			b1 = A[1][0] - centeriod[0]; b2 = A[1][1] - centeriod[1]; b3 = A[1][2] - centeriod[2];
			c1 = A[2][0] - centeriod[0]; c2 = A[2][1] - centeriod[1]; c3 = A[2][2] - centeriod[2];
		}
		else {
			a1 = A[0][0]; a2 = A[0][1]; a3 = A[0][2];
			b1 = A[1][0]; b2 = A[1][1]; b3 = A[1][2];
			c1 = A[2][0]; c2 = A[2][1]; c3 = A[2][2];
		}
		det = a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) + a3 * (b1 * c2 - b2 * c1);

		double term1 = det * factorial(p) * factorial(q) * factorial(r) / factorial(s + 3);
		double term2 = 0;
		vector<vector<double>> K;
		for (unsigned int i = 0; i < s_out.size(); i++) {
			K = s_out[i];
			double numerator = factorial(K[0][0] + K[1][0] + K[2][0]) * factorial(K[0][1] + K[1][1] + K[2][1]) * factorial(K[0][2] + K[1][2] + K[2][2]);
			double denominator = factorial(K[0][0]) * factorial(K[0][1]) * factorial(K[0][2]) * factorial(K[1][0]) * factorial(K[1][1]) * factorial(K[1][2]) * factorial(K[2][0]) * factorial(K[2][1]) * factorial(K[2][2]);
			double exterm_right_term = pow(a1, K[0][0]) * pow(b1, K[0][1]) * pow(c1, K[0][2]) * pow(a2, K[1][0]) * pow(b2, K[1][1]) * pow(c2, K[1][2]) * pow(a3, K[2][0]) * pow(b3, K[2][1]) * pow(c3, K[2][2]);
			term2 = term2 + ((numerator / denominator) * exterm_right_term);
		}
		moment = moment + (term1 * term2);
		}
	if (isScaled) {
		moment = moment / (pow(volume, (1 + ((s) / 3.0))));
	}

	return moment;
}

vector<pair<string, double>> geometricMoments::getMomentVector(double order, bool isCentered, bool isScaled) {
	vector<double> temp_1;
	vector<vector<double>> pqr;

	// getting all the moments of a specific order
	for (int a = 0; a <= 10; a++) {
		for (int b = 0; b <= 10; b++) {
			double t = order - a - b;
			vector<double> temp;
			if ((0 <= t) && (t < 10)) {
				temp_1.clear();
				temp_1.push_back(a);
				temp_1.push_back(b);
				temp_1.push_back(order - a - b);
				pqr.push_back(temp_1);
			}
		}
	}

	/// remember to remove these lines 
		// Sort the indices in descending order (required for proper erasing)
		//vector<int> indicesToErase = { 23,21,19,16,14,12,8,5,3,1 };
		////sort(indicesToErase.rbegin(), indicesToErase.rend());

		//// Erase the elements in reverse order to avoid invalidating subsequent indices
		//for (int index : indicesToErase) {
		//	pqr.erase(pqr.begin() + index);
		//}
	//--------------------------

	vector<pair<string, double>> momentVector;
	stringstream momenName;

	for (int i = 0; i < pqr.size(); i++) {
		// geting moment names as string
		momenName.str(string());
		momenName << "m_";
		for (int j = 0; j < pqr[i].size(); j++) {
			momenName << pqr[i][j];
		}
		// getting the moments 
		momentVector.push_back(make_pair(momenName.str(), getMoment(pqr[i][0], pqr[i][1], pqr[i][2], isCentered, isScaled)));
	}
	return momentVector;
}


void geometricMoments::getVolumeProperties(vector<double> &centeriod, double &volume) {
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, det;
	double m_000 = 0, m_010 = 0, m_100 = 0, m_001 = 0;
	for (size_t itri = 0; itri < m_mesh.num_tris(); ++itri) {
		vector<vector<double>> A;
		const float* c;
		vector<double> A_temp;

		c = m_mesh.tri_corner_coords(itri, 0);
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		c = m_mesh.tri_corner_coords(itri, 1);
		A_temp.clear();
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		c = m_mesh.tri_corner_coords(itri, 2);
		A_temp.clear();
		A_temp.push_back(c[0]); A_temp.push_back(c[1]); A_temp.push_back(c[2]);
		A.push_back(A_temp);

		a1 = A[0][0]; a2 = A[0][1]; a3 = A[0][2];
		b1 = A[1][0]; b2 = A[1][1]; b3 = A[1][2];
		c1 = A[2][0]; c2 = A[2][1]; c3 = A[2][2];
		det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

		m_000 = m_000 + (1.0 / 6.0 * det);
		m_100 = m_100 + (1.0 / 24.0 * det * (a1 + b1 + c1));
		m_010 = m_010 + (1.0 / 24.0 * det * (a2 + b2 + c2));
		m_001 = m_001 + (1.0 / 24.0 * det * (a3 + b3 + c3));
	}
	volume = m_000;
	centeriod.push_back(m_100 / m_000);
	centeriod.push_back(m_010 / m_000);
	centeriod.push_back(m_001 / m_000);
}

vector<vector<vector<double>>>	 geometricMoments::momentMatrix(double p, double q, double r) {
	vector<double> temp_1;
	vector<vector<double>> k1, k2, k3, temp_2;
	vector<vector<vector<double>>> K;

	for (int a = 0; a <= 10; a++) {
		for (int b = 0; b <= 10; b++) {
			double t = p - a - b;
			vector<double> temp;
			if ((0 <= t) && (t < 10)) {
				temp_1.clear();
				temp_1.push_back(a);
				temp_1.push_back(b);
				temp_1.push_back(p - a - b);
				k1.push_back(temp_1);
			}
		}
	}

	for (int a = 0; a <= 10; a++) {
		for (int b = 0; b <= 10; b++) {
			double t = q - a - b;
			if ((0 <= t) && (t < 10)) {
				temp_1.clear();
				temp_1.push_back(a);
				temp_1.push_back(b);
				temp_1.push_back(q - a - b);
				k2.push_back(temp_1);
			}
		}
	}

	for (int a = 0; a <= 10; a++) {
		for (int b = 0; b <= 10; b++) {
			double t = r - a - b;
			if ((0 <= t) && (t < 10)) {
				temp_1.clear();
				temp_1.push_back(a);
				temp_1.push_back(b);
				temp_1.push_back(r - a - b);
				k3.push_back(temp_1);
			}
		}
	}

	for (unsigned int i = 0; i < k1.size(); i++) {
		for (unsigned int j = 0; j < k2.size(); j++) {
			for (unsigned int k = 0; k < k3.size(); k++) {
				temp_2.clear();
				temp_2.push_back(k1[i]);
				temp_2.push_back(k2[j]);
				temp_2.push_back(k3[k]);
				K.push_back(temp_2);
			}
		}
	}
	return K;
}


double geometricMoments::factorial(double n) {
	if (n > 1)
		return n * factorial(n - 1);
	else
		return 1;
}