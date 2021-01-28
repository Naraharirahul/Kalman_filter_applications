#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0.02, 0,
      0, 0, 0, 0, 0.02;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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
  n_x_ = 5;
  lambda_ = 3 - n_x_;
//   cout << lambda_ << endl;
  n_aug_ = n_x_ + 2;
  double n_sig_ = 2 * n_aug_ + 1; 
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  time_us_ = 0.0;
  
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));  weights_(0) = lambda_ / (lambda_ + n_aug_) ;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //.. Initialize measurement noise covairance matix
  // int n_z_ = 3;
  // MatrixXd R = MatrixXd(n_z_, n_z_);  // For the radar

  // R << std_radr_ * std_radr_, 0, 0,
  //      0, std_radphi_*std_radphi_,0,
  //      0, 0, std_radrd_*std_radr_;
  
  // MatrixXd R_;   // For the Lidar
  // R_.setZero(2, 2);
  // R_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  
   if (!is_initialized_)
    {
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

        // P_ << std_radr_* std_radr_, 0, 0, 0, 0,
        //           0, std_radr_ * std_radr_, 0, 0, 0,
        //           0, 0, std_radrd_ * std_radrd_, 0, 0,
        //           0, 0, 0, std_radphi_ * std_radphi_, 0,
        //           0, 0, 0, 0, std_radphi_ * std_radphi_;  
    }
  else if (meas_package.sensor_type_ ==  MeasurementPackage::LASER){
    double x = meas_package.raw_measurements_[0];
    double y = meas_package.raw_measurements_[1];

    x_ << x, y, 0, 0, 0;

    // P_ << std_laspx_, 0, 0, 0, 0,
    //       0, std_laspy_, 0, 0, 0,
    //       0, 0, 1, 0, 0,
    //       0, 0, 0, 1, 0,
    //       0, 0, 0, 0, 1;  
  }
  is_initialized_ = true;
  return;

  }
  double time_diff = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(time_diff);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  VectorXd x_aug = VectorXd(n_aug_);

  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // x_aug.fill(0.0);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  Xsig_aug.col(0) = x_aug;

  MatrixXd L = P_aug.llt().matrixL();

  for(int i = 0; i <  n_aug_ ; i++){
//     cout<< i << "  Value of i" << endl;
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  for(int i = 0; i < 2 * n_aug_ + 1 ; i++){
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_p, py_p,v_p,yaw_p,yawd_p;

    if(fabs(yawd) > 0.001){
      px_p = p_x + v/yawd * (sin (yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yaw*delta_t));
    }
    else{
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    v_p = v;
    yaw_p = yaw + yawd*delta_t;
    yawd_p = yawd;

    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  x_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){

    x_ = x_ +  weights_(i) * Xsig_pred_.col(i);
  }
  P_.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++){

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3) > M_PI)
    {
      x_diff(3) -= 2.*M_PI;
    }

    while(x_diff(3) < -M_PI)
    {
      x_diff(3) += 2.*M_PI;
    }

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

  }


  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z_ = 2;
  MatrixXd R_;   // For the Lidar
  R_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

  MatrixXd Zsig_ = MatrixXd(n_z_, 2*n_aug_+1);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    Zsig_(0, i) = Xsig_pred_(0, i);
    Zsig_(1, i) = Xsig_pred_(1, i);
  }

  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //.. calculate Kalman gain K
  MatrixXd K = Tc * S.inverse();

  VectorXd z_ = meas_package.raw_measurements_;
  VectorXd z_diff = z_ - z_pred_;
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose(); 

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z_ = 3;
  VectorXd z_ = meas_package.raw_measurements_;
  MatrixXd z_sig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);

  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    z_sig(0,i) = sqrt(p_x * p_x + p_y*p_y);
    z_sig(1,i) = atan2(p_y,p_x);
    z_sig(2,i) = (p_x * v1 + p_y*v2) / (sqrt(p_x*p_x + p_y*p_y));
  }

  z_pred.fill(0.0);

  for(int i = 0; i < 2* n_aug_ + 1; i++){
    z_pred = z_pred + weights_(i) * z_sig.col(i); 
    }

  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  
  for(int i =0 ; i < 2* n_aug_ + 1; i++){

    VectorXd z_diff = z_sig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_, n_z_);

  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_*std_radphi_,0,
       0, 0, std_radrd_*std_radr_;

  S = S+R;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  for(int i = 0; i < 2* n_aug_ +1 ; i++){


    VectorXd z_diff = z_sig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;


    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff =  z_ - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}