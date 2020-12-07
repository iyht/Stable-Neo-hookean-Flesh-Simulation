#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    Eigen::Matrix1212d M;
    mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

    Eigen::Vector12d qdot_el;
    qdot_el << qdot.segment(3*element(0),3),
               qdot.segment(3*element(1),3),
               qdot.segment(3*element(2),3),
               qdot.segment(3*element(3),3);
    T = 0.5*qdot_el.transpose()*M*qdot_el;

}