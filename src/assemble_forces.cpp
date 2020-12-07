#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D, bool stable) {


    int r = q.rows();
    f.resize(r);
    f.setZero();

    for(int i = 0; i < T.rows(); i++)
    {
        Eigen::Vector3i index_q0, index_q1, index_q2, index_q3; // store the index to the vector q to reach [q0_x, q0_y, q0_z]
        index_q0 << T(i, 0)*3, T(i, 0)*3+1, T(i, 0)*3+2;
        index_q1 << T(i, 1)*3, T(i, 1)*3+1, T(i, 1)*3+2;
        index_q2 << T(i, 2)*3, T(i, 2)*3+1, T(i, 2)*3+2;
        index_q3 << T(i, 3)*3, T(i, 3)*3+1, T(i, 3)*3+2;


        Eigen::Vector12d dV_i;
        dV_linear_tetrahedron_dq(dV_i, q, V, T.row(i), v0(i), C, D, stable);

        f(index_q0(0)) -= dV_i(0);
        f(index_q0(1)) -= dV_i(1);
        f(index_q0(2)) -= dV_i(2);

        f(index_q1(0)) -= dV_i(3);
        f(index_q1(1)) -= dV_i(4);
        f(index_q1(2)) -= dV_i(5);

        f(index_q2(0)) -= dV_i(6);
        f(index_q2(1)) -= dV_i(7);
        f(index_q2(2)) -= dV_i(8);

        f(index_q3(0)) -= dV_i(9);
        f(index_q3(1)) -= dV_i(10);
        f(index_q3(2)) -= dV_i(11);
    }



};