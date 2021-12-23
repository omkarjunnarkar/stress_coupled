

using namespace std;
using namespace Eigen;

double sumofsquaresofcoeff(MatrixXd P);
MatrixXd function_y(MatrixXd para);
double rho_function(MatrixXd para, MatrixXd y, MatrixXd y_plus_h, MatrixXd h, MatrixXd Jacobian, MatrixXd y_measured, double lambda);
MatrixXd jacobian_function(MatrixXd para, MatrixXd y, MatrixXd initial_deflection);
MatrixXd levenberg_fit(MatrixXd y_measured, MatrixXd initial_deflection, MatrixXd para, double lambda);

