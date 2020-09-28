#include <iostream>
// #include "Dense"
#include "ukf.h"
#include <eigen3/Eigen/Core>
using Eigen::MatrixXd;

int main() {

  // Create a UKF instance
  UKF ukf;

  /**
   * Programming assignment calls
   */
  MatrixXd Xsig = MatrixXd(7, 15);
  ukf.GenerateSigmaPoints(&Xsig);

  // print result
  std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  return 0;
}