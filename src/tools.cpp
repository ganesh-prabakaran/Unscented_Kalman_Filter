#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
  std::cout << "CalculateRMSE" << std::endl;
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
      * Calculate the RMSE here.
  */
  //std::cout << "CalculateRMSE" << std::endl;
  VectorXd  rmse(4);
  rmse << 0,0,0,0;
  if (estimations.size() != ground_truth.size() || estimations.size() == 0){
    return rmse;
  }
  //std::cout << "Estimations Size" << estimations.size() << endl;
  for(unsigned int i = 0; i < estimations.size(); ++i ){
    //std::cout << "EST" << i << estimations[i] << endl;
    //std::cout << "GT" << i << ground_truth[i] << endl;
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse +=  residual;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  //std::cout << "CalculateRMSE completed" << std::endl;
  //std::cout << "rmse :" ;
  //std::cout << rmse << std::endl;
  return rmse;
}
