#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
    const vector<VectorXd> &ground_truth) {
  int e_size = estimations.size();

  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;
  if (e_size < 1 || e_size != ground_truth.size()) return rmse;

  VectorXd summ(4);
  summ << 0, 0, 0, 0;
  for (int i=0; i<e_size; i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / e_size;
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if(px==0 && py==0){
    return Hj;
  }

  float px2 = px*px;
  float py2 = py*py;
  float pnorm2 = px2 + py2;
  float pnorm = std::sqrt(pnorm2);

  //compute the Jacobian matrix
  Hj(0,0) = px/pnorm;
  Hj(0,1) = py/pnorm;

  Hj(1,0) = -1 * py/pnorm2;
  Hj(1,1) = px/pnorm2;

  Hj(2,0) = py*(vx*py - vy*px)/(pnorm2*pnorm);
  Hj(2,1) = px*(vy*px - vx*py)/(pnorm2*pnorm);
  Hj(2,2) = px/pnorm;
  Hj(2,3) = py/pnorm;
  return Hj;
}
