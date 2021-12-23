#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include<math.h>
#include<vector>
#include"src/rapidcsv.h"
#include"feamain.h"
#include"LevenbergMarquardtFit.h"

using namespace std;
using namespace Eigen;
using namespace rapidcsv;

void main() {
	
	cout << "Starting the LMA Main program" << endl;
	
	Document dy("stress.csv");
	vector<double> ycol = dy.GetColumn<double>(0);
	int data_size = ycol.size();
	MatrixXd y_measured(ycol.size(), 1);

	for (int u = 0; u < data_size ; u++) {
		y_measured(u, 0) = ycol[u];
	}
	
	MatrixXd parameters{
		{90007} 
	};

	MatrixXd init{
		{0.1}
	};

	parameters = parameters.reshaped(parameters.cols(), 1);
	init = init.reshaped(parameters.cols(), 1);

	cout << "All Data Read, Calling LMA for Parameters" << endl;
	MatrixXd final_parameters = levenberg_fit(y_measured, init, parameters, 1e4);
	cout << final_parameters << endl;
}