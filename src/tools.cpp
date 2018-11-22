#include <iostream>
#include "tools.h"


#define PI 3.14159265

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
  //pre-compute a set of terms to avoid repeated calculation

  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero

  if(fabs(c1) < 0.0001){
    return Hj;
  }

  //compute the Jacobian matrix
  Hj <<  (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::Polar2Cartesian(const VectorXd& x_state){
  MatrixXd p2c(4, 3);
  float s = sin(x_state(1));
  float c = cos(x_state(1));
  p2c << c, 0, 0,
         s, 0, 0,
         0, 0, c,
         0, 0, s;
  return p2c * x_state;
}

VectorXd Tools::Cartesian2Polar(const VectorXd& x_state){
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float norm2 = px*px + py*py;
  float norm = std::sqrt(norm2);
  float projection = px*vx + py*vy;

  float theta = std::atan2(py, px);

  VectorXd polar(3);
  polar(0) = norm;
  polar(1) = theta;
  polar(2) = projection/norm;

  return polar;
}

float Tools::NormalizeAngle(float angle){
  float newAngle = angle;
  if(angle > PI) {
    newAngle -= 2*PI;
  } else if(angle < -1*PI) {
    newAngle += 2*PI;
  }
  return newAngle;
}
