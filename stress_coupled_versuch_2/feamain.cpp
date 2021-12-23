#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"
#include"feamain.h"

using namespace std;
using namespace Eigen;

MatrixXd feamain(MatrixXd parameters) {

	//ofstream mydisplacefile("displacement.csv");
	//ofstream mystressfile("stress.csv");
	//ofstream mystrainfile("strain.csv");

	double E = parameters(0, 0);
	//double sigma_0 = parameters(1, 0)*1e2;
	//double eta = parameters(1, 0);

	double F_max = 3900;
	double A_1 = 6;
	double A_2 = 12;
	double L_1 = 20;
	double L_2 = 40;
	//double E = 100000;
	int steps_num = 10000;

	MatrixXd output_stress(steps_num+2, 1);

	double steps = steps_num;
	double t_tot = 0.5;
	double sigma_0 = 200;
	int elements = 2;
	double force = 0;
	int gauss_weights = 2;
	double eta = 1;
	MatrixXd Areas, Lengths;

	steps = t_tot / double(steps);

	tie(Areas, Lengths) = get_areas_lengths(elements, A_1, A_2, L_1, L_2);
	//cout << "Areas=" << Areas << endl;
	//cout << "Lengths=" << Lengths << endl;

	MatrixXd F_ext = MatrixXd::Zero(elements + 1,1);
	MatrixXd E_mat= E*MatrixXd::Ones(elements , 1);
	MatrixXd eps_p = MatrixXd::Zero(elements , 1);
	MatrixXd eps = MatrixXd::Zero(elements , 1);
	MatrixXd u = MatrixXd::Zero(elements + 1, 1);
	MatrixXd K_t, F_int, strain, stress;
	
	int count = 0;

	while (force < F_max) {
		if (count != 0) { force = force + (steps / t_tot) * F_max; }
		
		count++;
		
		//cout << "Force = " << force << endl;

		F_ext(int(elements / 2), 0) = force;

		for (int NR = 0; NR < 6; NR++) {

			tie(K_t, F_int, eps_p, strain, stress) = material_routine(elements, gauss_weights, Lengths, Areas, E_mat, u, steps, eta, eps_p, sigma_0, E);
			//cout << "K_t= \n" << K_t << endl;
			MatrixXd K_t_red = K_t.block(1, 1, K_t.cols() - 2, K_t.rows() - 2);
			//cout << "K_tred= \n" << K_t_red << endl;
			MatrixXd K_t_red_inv = K_t_red.inverse();

			MatrixXd R = F_int - F_ext;
			//cout << "R=" << R << endl;
			MatrixXd R_red(elements - 1, 1),u_red(elements-1,1);
			for (int c = 1; c < elements; c++) { R_red(c - 1, 0) = R(c, 0); }
			//cout << "R red=" << R_red << endl;
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }
			//cout << "u_red=" << u_red << endl;
			MatrixXd del_u_red = K_t_red_inv * R_red;

			//cout << "del_u_red=" << del_u_red << endl;

			u_red = u_red - del_u_red;
			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }
			//cout << "u=" << u;

			//if (del_u_red.norm() < 1e-4) { break; }
		}
		
		//mydisplacefile << u((elements / 2) , 0) << endl;
		//mystressfile << stress((elements / 2)-1, 0) << endl;
		//mystrainfile << strain((elements / 2) - 1, 0) << endl;
		
		output_stress(count - 1, 0) = stress((elements / 2) - 1, 0);
	}

	//mydisplacefile.close();
	//mystressfile.close();
	//mystrainfile.close();

	return output_stress;
}