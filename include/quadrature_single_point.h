#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {

    // q element
    Eigen::Vector12d q_el;
    q_el << q.segment(3*element(0),3), q.segment(3*element(1),3), q.segment(3*element(2),3), q.segment(3*element(3),3);

    // Centroid of a tet
    Eigen::Vector3d X;
    X[0] = (q_el[0] + q_el[3] + q_el[6] + q_el[9])/4;  // x
    X[1] = (q_el[1] + q_el[4] + q_el[7] + q_el[10])/4; // y
    X[2] = (q_el[2] + q_el[5] + q_el[8] + q_el[11])/4; // z

    //double out;
    integrand(integrated, q, element, X);

    integrated = integrated*volume;

}

