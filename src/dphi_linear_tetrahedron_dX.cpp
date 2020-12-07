#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Vector3d X_0, X_1, X_2, X_3;
    X_0 = V.row(element(0));
    X_1 = V.row(element(1));
    X_2 = V.row(element(2));
    X_3 = V.row(element(3));

    Eigen::Matrix3d T, T_inv;
    T.setZero();
    T.block(0, 0, 3, 1) = X_1- X_0;
    T.block(0, 1, 3, 1) = X_2- X_0;
    T.block(0, 2, 3, 1) = X_3- X_0;

    T_inv = T.inverse();

    dphi.setZero();
    dphi.block(0, 0, 1, 3) = -T_inv.colwise().sum();
    dphi.block(1, 0, 3, 3) = T_inv;

}