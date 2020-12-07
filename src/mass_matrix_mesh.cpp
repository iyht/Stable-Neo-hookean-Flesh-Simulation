#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    int r = qdot.rows(); // r = 3n
    M.setZero();
    M.resize(r, r);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(r*r);

    for(int i = 0; i < T.rows(); i++)
    {

        Eigen::Matrix1212d M_i;
        mass_matrix_linear_tetrahedron(M_i, qdot, T.row(i), density, v0(i));


        Eigen::Matrix3d M_i00 = M_i.block(0, 0, 3, 3);
        Eigen::Matrix3d M_i01 = M_i.block(0, 3, 3, 3);
        Eigen::Matrix3d M_i02 = M_i.block(0, 6, 3, 3);
        Eigen::Matrix3d M_i03 = M_i.block(0, 9, 3, 3);

        Eigen::Matrix3d M_i10 = M_i.block(3, 0, 3, 3);
        Eigen::Matrix3d M_i11 = M_i.block(3, 3, 3, 3);
        Eigen::Matrix3d M_i12 = M_i.block(3, 6, 3, 3);
        Eigen::Matrix3d M_i13 = M_i.block(3, 9, 3, 3);

        Eigen::Matrix3d M_i20 = M_i.block(6, 0, 3, 3);
        Eigen::Matrix3d M_i21 = M_i.block(6, 3, 3, 3);
        Eigen::Matrix3d M_i22 = M_i.block(6, 6, 3, 3);
        Eigen::Matrix3d M_i23 = M_i.block(6, 9, 3, 3);

        Eigen::Matrix3d M_i30 = M_i.block(9, 0, 3, 3);
        Eigen::Matrix3d M_i31 = M_i.block(9, 3, 3, 3);
        Eigen::Matrix3d M_i32 = M_i.block(9, 6, 3, 3);
        Eigen::Matrix3d M_i33 = M_i.block(9, 9, 3, 3);

        for(int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {

                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 0)+jj, M_i00.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 1)+jj, M_i01.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 2)+jj, M_i02.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 0)+ii, 3*T(i, 3)+jj, M_i03.coeff(ii, jj)));


                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 0)+jj, M_i10.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 1)+jj, M_i11.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 2)+jj, M_i12.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 1)+ii, 3*T(i, 3)+jj, M_i13.coeff(ii, jj)));

                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 0)+jj, M_i20.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 1)+jj, M_i21.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 2)+jj, M_i22.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 2)+ii, 3*T(i, 3)+jj, M_i23.coeff(ii, jj)));

                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 0)+jj, M_i30.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 1)+jj, M_i31.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 2)+jj, M_i32.coeff(ii, jj)));
                tripleList.push_back(Trip(3*T(i, 3)+ii, 3*T(i, 3)+jj, M_i33.coeff(ii, jj)));

            }
        }

    }
    M.setFromTriplets(tripleList.begin(), tripleList.end());
}
