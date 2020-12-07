#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D, bool stable) {

    int r = qdot.rows(); // r = 3n
    K.setZero();
    K.resize(r, r);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(r*r);

    for(int i = 0; i < T.rows(); i++)
    {

        Eigen::Matrix1212d H_i;
        d2V_linear_tetrahedron_dq2(H_i, q, V, T.row(i), v0(i), C, D, stable);


        Eigen::Matrix3d H_i00 = H_i.block(0, 0, 3, 3);
        Eigen::Matrix3d H_i01 = H_i.block(0, 3, 3, 3);
        Eigen::Matrix3d H_i02 = H_i.block(0, 6, 3, 3);
        Eigen::Matrix3d H_i03 = H_i.block(0, 9, 3, 3);

        Eigen::Matrix3d H_i10 = H_i.block(3, 0, 3, 3);
        Eigen::Matrix3d H_i11 = H_i.block(3, 3, 3, 3);
        Eigen::Matrix3d H_i12 = H_i.block(3, 6, 3, 3);
        Eigen::Matrix3d H_i13 = H_i.block(3, 9, 3, 3);

        Eigen::Matrix3d H_i20 = H_i.block(6, 0, 3, 3);
        Eigen::Matrix3d H_i21 = H_i.block(6, 3, 3, 3);
        Eigen::Matrix3d H_i22 = H_i.block(6, 6, 3, 3);
        Eigen::Matrix3d H_i23 = H_i.block(6, 9, 3, 3);

        Eigen::Matrix3d H_i30 = H_i.block(9, 0, 3, 3);
        Eigen::Matrix3d H_i31 = H_i.block(9, 3, 3, 3);
        Eigen::Matrix3d H_i32 = H_i.block(9, 6, 3, 3);
        Eigen::Matrix3d H_i33 = H_i.block(9, 9, 3, 3);

        for(int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {

                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 0)+jj, -H_i00.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 1)+jj, -H_i01.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 2)+jj, -H_i02.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 3)+jj, -H_i03.coeff(ii, jj)));


                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 0)+jj, -H_i10.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 1)+jj, -H_i11.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 2)+jj, -H_i12.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 3)+jj, -H_i13.coeff(ii, jj)));

                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 0)+jj, -H_i20.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 1)+jj, -H_i21.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 2)+jj, -H_i22.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 3)+jj, -H_i23.coeff(ii, jj)));

                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 0)+jj, -H_i30.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 1)+jj, -H_i31.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 2)+jj, -H_i32.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 3)+jj, -H_i33.coeff(ii, jj)));

            }
        }

    }
    K.setFromTriplets(tripleList.begin(), tripleList.end());
}
