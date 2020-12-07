#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {


    int l = V_skin.rows();
    int n = V.rows();
    N.resize(l, n);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(l*n);

    for(int i=0; i < l; i ++)
    {
        Eigen::Vector3d x = V_skin.row(i);
        bool mapped = false;
        double min_dist = std::numeric_limits<int>::max();
        Eigen::Vector3d min_v;
        Eigen::Vector4d min_phi;
        Eigen::Vector4i min_element;

        // search all the vertices on tet that has the min_dist with current mesh vertex
        for(int j = 0; j < n; j ++)
        {
            Eigen::Vector4d phi;
            Eigen::RowVector4i element = T.row(j);
            phi_linear_tetrahedron(phi, V, element, x);
            for(int v = 0; v < 4; v++)
            {
                Eigen::Vector3d X_v;
                X_v = V.row(element(v));
                double dist = (X_v - x).norm();
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_phi = phi;
                    min_v = X_v;
                    min_element = element;
                }

            }
        }
        tripleList.push_back(Trip(i, min_element[0], min_phi[0]));
        tripleList.push_back(Trip(i, min_element[1], min_phi[1]));
        tripleList.push_back(Trip(i, min_element[2], min_phi[2]));
        tripleList.push_back(Trip(i, min_element[3], min_phi[3]));
    }
    N.setFromTriplets(tripleList.begin(), tripleList.end());

}