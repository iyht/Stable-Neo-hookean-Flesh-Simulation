#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {

    double p = 0.5; // scaling factor
    double c = 1e-8; // ensure sufficient decrease
    //solve the linear system H*d = -g to the d
    // A = H
    // b = -g
    for(int i = 0; i < maxSteps; i++)
    {
        g(tmp_g, x0);
        H(tmp_H, x0);
        if(tmp_g.norm() < c)
        {
            break;
        }
        //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cgSolver;
        //cgSolver.compute(tmp_H);
        //Eigen::VectorXd d;
        //d = cgSolver.solve(-tmp_g);

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(tmp_H);
        Eigen::VectorXd d;
        d = solver.solve(-tmp_g);

        //std::cout << "d\n" << d <<  std::endl;
        //std::cout << "x0\n" << x0 <<  std::endl;

        // line search
        double alpha = 1.0; // initial step length
        while(true)
        {
           double threshold_cost = f(x0) + c*d.transpose()*tmp_g;
           double new_cost = f(x0 + (alpha*d));
           if(new_cost <= threshold_cost || alpha < c)
           {
               //std::cout << "alpha\n" << alpha <<  std::endl;
               break;
           }
           alpha = p*alpha;
        }
        x0 = x0 + alpha*d;
    }

   return 0.0;
}
