#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .2;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
    Initialization - State, Process, Measurement Noise Covariance matrix
   */
  time_us_ = 0.0;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  Q_ = MatrixXd(2,2);
  Q_ << std_a_ * std_a_ , 0.0,
        0.0, std_yawdd_ * std_yawdd_;

  R_radar = MatrixXd(3,3);
  R_radar << std_radr_ * std_radr_ , 0.0, 0.0,
        0.0, std_radphi_ * std_radphi_, 0.0,
        0.0, 0.0, std_radrd_ * std_radrd_;

  R_lidar = MatrixXd(2,2);
  R_lidar << std_laspx_ * std_laspx_ , 0.0,
        0.0, std_laspy_ * std_laspy_;

  /** Initialize Weights **/
  weights_ = VectorXd( 2 * n_aug_ + 1);
  weights_[0]  = lambda_ / (lambda_ + n_aug_);
  for (unsigned int i = 1; i < 2 * n_aug_ + 1; i++){
      weights_[i] = 0.5 / (lambda_ + n_aug_);
  }

  /** Initialize State Covariance Matrix **/
  P_ << 0.15, 0,   0, 0, 0,
        0,    0.15,0, 0, 0,
        0,    0,   1, 0, 0,
        0,    0,   0, 1, 0,
        0,    0,   0, 0, 1;
  /** Define State prediction matrix for sigma points **/
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_){
     x_ << 1, 1, 1, 1, 0.1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      float p_x = meas_package.raw_measurements_[0];
      float p_y = meas_package.raw_measurements_[1];
      x_[0] = p_x;
      x_[1] = p_y;

      if (fabs(p_x) < 0.001){
      x_[0] = 0.001;
      }
      if (fabs(p_y) < 0.001){
      x_[1] = 0.001;
      }
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ){
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];
      float vx = ro_dot * cos(theta);
      float vy = ro_dot * sin(theta);
      float v  = sqrt(vx * vx + vy * vy);
      x_[0] = ro * cos(theta);
      x_[1] = ro * sin(theta);

    }

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;

  }

  long long delta_t_us = meas_package.timestamp_ - time_us_;
  double delta_t = static_cast<double>(delta_t_us) / static_cast<double>(1e6) ;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
    // Radar updates
    UpdateRadar(meas_package);

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER ){
    // Laser updates
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

    //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);

  //create sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, 2*n_x_+1);

  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q_;

  //calculate square root of P
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  //cout << "Updated col 0 of Xsig_aug" << endl;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

 //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  // Predicted state mean
  x_.fill(0.0);

  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i ++){
      x_ = x_ + weights_[i] * Xsig_pred_.col(i);
  }

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i ++){
      VectorXd x_diff =  Xsig_pred_.col(i)  - x_;

      //angle normalization
      x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));

      P_ = P_ + (weights_[i] * x_diff * x_diff.transpose());
  }
}
/**
 * Updates the state and the state covariance matrix using a laser/Radar measurement.
 * @param {MeasurementPackage} meas_package

 Calculate the NIS for both lidar and radar.
 */
void UKF::Update_x_P( MatrixXd Zsig, VectorXd z, MatrixXd R ){
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if (n_z_ == 3) {
            z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    }
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R;

  /////////
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  //calculate cross correlation matrix

  for ( unsigned int i = 0; i<2*n_aug_+1; i++){
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //State difference
      VectorXd x_diff = Xsig_pred_.col(i)  - x_;
      if (n_z_ == 3){
      //angle normalization
      z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
      }

      //angle normalization
      x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z_diff = z - z_pred;
  if (n_z_ == 3){ // Radar
    // Angle normalization
      z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
  }

  MatrixXd K = MatrixXd(n_x_, n_z_);
  //calculate Kalman gain K;
  K = Tc * S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();

  // Calculate NIS and write into file
 NIS_ = z_diff.transpose() * S.inverse() * z_diff;
 ofstream outfile;
 outfile.open ("outputNIS.txt", ios::out | ios::app );
 outfile << n_z_ << " " << NIS_ << endl;
 outfile.close();
}


/**
 * Generate updates for state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.


  */

  n_z_ = 2;
  double p_x = meas_package.raw_measurements_[0];
  double p_y = meas_package.raw_measurements_[1];

  VectorXd z = VectorXd(n_z_);

  z << p_x, p_y;

  //MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig = Xsig_pred_.block(0, 0, n_z_, 2 * n_aug_ + 1);

  Update_x_P(Zsig, z,  R_lidar);
}


/**
 * Generate updates for state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  */
  n_z_ = 3;
  double ro = meas_package.raw_measurements_[0];
  double theta = meas_package.raw_measurements_[1];
  double ro_dot = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(n_z_);
  z << ro,theta, ro_dot;
  //cout << "Retrieve radar measurements" << z << endl;

    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
  for (int i=0; i<2*n_aug_+1; i++) {  //2n+1 sigma points

    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
  }
  Update_x_P(Zsig, z, R_radar);

}
