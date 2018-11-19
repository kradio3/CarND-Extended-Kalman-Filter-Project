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
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  return MatrixXd(7,7);
}
