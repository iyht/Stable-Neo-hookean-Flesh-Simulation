#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {


    int l = indices.size();
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;
    tripleList.reserve(q_size*q_size);

    int j=0;
    for(int i=0; i<q_size; i+=3) {

        if(j < indices.size() && indices[j] == i/3)
        {
            j += 1;
        }
        else
        {
            tripleList.push_back(T(i-3*j +0 ,i+0, 1.));
            tripleList.push_back(T(i-3*j +1 ,i+1, 1.));
            tripleList.push_back(T(i-3*j +2 ,i+2, 1.));
        }

    }
    P.resize(q_size - (3*l), q_size);
    P.setFromTriplets(tripleList.begin(), tripleList.end());


}