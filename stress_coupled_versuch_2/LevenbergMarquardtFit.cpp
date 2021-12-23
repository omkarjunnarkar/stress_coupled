
#include<iostream>
#include<Eigen/Dense>
#include"feamain.h"
#include"LevenbergMarquardtFit.h"

using namespace std;
using namespace Eigen;

double sumofsquaresofcoeff(MatrixXd P) {
	double s = 0;

	for (int i = 0; i < P.rows(); i++) {
		for (int j = 0; j < P.cols(); j++) {
			s += P(i, j) * P(i, j);
		}
	}
	return s;
};

/*To compute the value of function y=f(x) ; arguments : Parameters, x */
MatrixXd function_y(MatrixXd para) {
	MatrixXd y = feamain(para);
	return y;
};

/*To compute the deciding factor 'rho' for increase or decrease of damping factor Lambda*/
double rho_function(MatrixXd para, MatrixXd y, MatrixXd y_plus_h, MatrixXd h, MatrixXd Jacobian, MatrixXd y_measured, double lambda) {
	//cout << "para=\n" << para << endl;
	//cout << "h=\n" << h << endl;
	//cout << "J=\n" << Jacobian << endl;
	//cout << "lambda=\n" << lambda << endl;

	MatrixXd value_arr = h.transpose() * (lambda * h + Jacobian.transpose() * (y_measured - y));
	double value = value_arr(0, 0);
	return (sumofsquaresofcoeff(y_measured - y) - sumofsquaresofcoeff(y_measured - y_plus_h)) / (para.rows() * value);
};

/*To compute the Jacobian Matrix by finite differences*/
MatrixXd jacobian_function(MatrixXd para, MatrixXd y, MatrixXd initial_deflection) {

	MatrixXd Jacobian_Matrix(y.rows(), para.rows());
	//MatrixXd y_para = function_y(para);
	MatrixXd y_deflected(y.rows(), 1);

	for (int i = 0; i < para.rows(); i++) {

		para(i, 0) = para(i, 0) + initial_deflection(i, 0);		/*Changing the parameters one by one */

		y_deflected = function_y(para);	/*Computing the deflected function arrray */
		for (int j = 0; j < y.rows(); j++) {
			Jacobian_Matrix(j, i) = (y_deflected(j, 0) - y(j, 0)) / initial_deflection(i, 0);
		}
		para(i, 0) = para(i, 0) - initial_deflection(i, 0);		/*Bringing back the parametes to original value*/
	}
	return Jacobian_Matrix;
};

/*Levenberg-Marquardt Algorithm ; Arguements: Actual function value, x ,Initial deflection of parameters , Initial guess of parameters , Damping factor Lambda */
MatrixXd levenberg_fit(MatrixXd y_measured, MatrixXd initial_deflection, MatrixXd para, double lambda) {

	cout << "Entered levenberg_fit" << endl;
	MatrixXd IdentityMat = MatrixXd::Identity(para.rows(), para.rows());
	int count = 0;
	double v = 2.0;
	double chi2p = sumofsquaresofcoeff(y_measured - function_y(para)) / para.rows(); /*Residual for initial guess of parameters*/
	cout << "Before do-while\n" << endl;

	do {
		count++;
		cout << "\nCount = " << count << endl;
		cout << "Inside do-while" << endl;
		//cout << "Parameters=\n" << para << endl;

		MatrixXd y = function_y(para);	/*Function Array*/
		cout << "Computed y: " << endl;
		double chi2p_plus_h;		/*Residual for next iteration*/
		MatrixXd J = jacobian_function(para, y, initial_deflection);	/*Jacobian Matrix*/
		//cout << "JtJ" << J.transpose() * J << endl;
		cout << "LAMBDA = " << lambda << endl;
		MatrixXd intermediate = J.transpose() * J + lambda * IdentityMat;
		//cout << "Intermediate " << intermediate << endl;
		MatrixXd h = intermediate.inverse() * (J.transpose() * (y_measured - y));	/*Computed Change in parameters*/
		cout << "Computed 'h'" << endl << h << endl;
		MatrixXd y_plus_h = function_y(para + h);
		//cout << "h =\n" << h << endl;
		cout << "Computed y_plus_h" << endl;
		chi2p = sumofsquaresofcoeff(y_measured - y) / para.rows();
		chi2p_plus_h = sumofsquaresofcoeff(y_measured - y_plus_h) / para.rows();
		cout << "Computed both chi squares" << endl;
		double rho = rho_function(para, y, y_plus_h, h, J, y_measured, lambda);	/*Computing the deciding factor for Lambda change*/
		cout << "Computed Rho" << endl;
		cout << "chi2p = " << chi2p << ", rho = " << rho << endl;

		if (rho > 1e-1) {
			para = para + h;
			v = 2.0;
			if (((1 - pow(((2 * rho)-1), 3)) < 1 / 3) && ((1 - pow(((2 * rho) - 1), 3))==0)) { lambda = lambda * (1 / 3); }
			else {
				lambda = lambda * (1 - pow(((2 * rho) - 1), 3));
			}

			cout << "Parameters updated" << endl;
			cout << "NOW LAMBDA = " << lambda << endl;
		}
		else {
			
			lambda = lambda * v;
			v = 2 * v;
			cout << "Parameters not updated" << endl;
			cout << "NOW LAMBDA = " << lambda << endl;
		}

		cout << "Iteration complete; "<<para << endl;
	} while ((chi2p * para.rows() / (para.rows() + y_measured.rows() - 1)) > 1e-2);
	//while ((chi2p*para.rows()/(para.rows() + x.rows() - 1)) > 1e-2);

	return para;
};