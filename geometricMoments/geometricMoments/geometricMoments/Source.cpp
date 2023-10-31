#include <iostream>
#include "geometricMoments.h"
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <math.h>
#include <cmath>
#include <fstream> 
#include <algorithm>
#include <atlstr.h>
#include "stl_reader.h"

using namespace std;

int main() {
	char filename[256];
	string filePath_read = "H:/Project with Prof Kaklis/01 BAR/igs_inDataset/stls/";
	string filePath_save = "H:/Project with Prof Kaklis/01 BAR/igs_inDataset/operators/";
	string init1 = ".stl", init3 = ".csv";
	int index2;

	// because all the input files have a different names we will those first. 
	// Read from the text file
	ifstream MyReadFile("H:/Project with Prof Kaklis/01 BAR/fileNamesForRhino.csv");
	vector<string> fileNamesForRhino; int i = 0;
	// Use a while loop together with the getline() function to read the file line by line
	string myText;
	while (getline(MyReadFile, myText)) {
		// Output the text from the file
		fileNamesForRhino.push_back(myText);
		i += 0;
	}
	vector<pair<string, double>> momentVector, momentVector_temp;
	for ( index2 = 0; index2 < fileNamesForRhino.size(); index2++) {
		cout << "-------------- Design: " << index2 << "\n";
		clock_t begin = clock();
			sprintf_s(filename, "%s%d", "_", index2);

			// Appending the string.
			string inputFile = filePath_read + fileNamesForRhino[i] + init1;

			//Output File Name
			sprintf_s(filename, "%s%d", "_", index2);
			string fileMoments = filePath_save + fileNamesForRhino[i] + init3;
			
			geometricMoments GM(inputFile);

			momentVector = GM.getMomentVector(0, 1, 0);
			//-- If you are calculationg multiple orders then uncomment it 
			for (int i= 1; i <= 4; i++) {
					momentVector_temp = GM.getMomentVector(i, 1, 0);
					momentVector.insert(momentVector.end(), momentVector_temp.begin(), momentVector_temp.end());
			}

			//saving to the file
			std::ofstream myfile;
			myfile.open(fileMoments);
			myfile << setprecision(16);
			for (unsigned int i = 0; i < momentVector.size(); i++) {
				myfile << momentVector[i].first << "," << momentVector[i].second << '\n';
			}
			myfile.close();

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "time: " << elapsed_secs << '\n';

	}




 	system("pause");
	return 0;
}
