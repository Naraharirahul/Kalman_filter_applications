#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  lambda_ = 3 - n_x_;
  n_aug_ = n_x_ + 2;
  double n_sig_ = 2 * n_aug_ + 1; 
  // std::cout << MeasurementPackage::raw_measurements_ << std::endl; 


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if(meas_package.sensor_type_ ==  MeasurementPackage::RADAR){
    
    // std::cou << meas_package.raw_measurements_ << std::endl; 
   double rho = meas_package.raw_measurements_[0]; 
   double phi = meas_package.raw_measurements_[1];
   double rho_dot = meas_package.raw_measurements_[2];

   double x = rho * cos(phi);
   double y = rho * sin(phi);
   double vx = rho_dot * cos(phi);
   double vy = rho_dot * sin(phi);
   double v = sqrt(vx*vx + vy*vy);

  x_ << x, y, v, rho, rho_dot ;

  P_ << std_radr_* std_radr_, 0, 0, 0, 0,
            0, std_radr_ * std_radr_, 0, 0, 0,
            0, 0, std_radrd_ * std_radrd_, 0, 0,
            0, 0, 0, std_radphi_ * std_radphi_, 0,
            0, 0, 0, 0, std_radphi_ * std_radphi_;  
  }
  else if (meas_package.sensor_type_ ==  MeasurementPackage::LASER){
    double x = meas_package.raw_measurements_[0];
    double y = meas_package.raw_measurements_[1];

    x_ << x, y, 0, 0, 0;

    P_ << std_laspx_, 0, 0, 0, 0,
          0, std_laspy_, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1; 
  }
    
  }
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */


void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}