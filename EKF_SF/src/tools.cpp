#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj = MatrixXd(3,4);

  double c1 = x_state(0);
  double c2 = x_state(1);
  double c3 = c1*c1 + c2*c2;
  double c4 = sqrt(c3);
  double c5 = x_state(2) * x_state(1) - x_state(3) * x_state(0);
  double c6 = -c5;

  Hj << c1/c4, c2/c4, 0, 0,
        -c2/c3, c1/c3, 0, 0,
        c2*c5/(c3*c4), c1*c6/c3*c4, c1/c4, c2/c4;

   return Hj;
}
