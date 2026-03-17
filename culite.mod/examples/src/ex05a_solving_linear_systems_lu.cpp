/**
 * @example ex05a_solving_linear_systems_lu.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/algebra.hpp>
#include <culite/dense.hpp>
#include <culite/linsol.hpp>

int main()
{
    /*
     * Create a random dense objects on host
     */
    const cla3p::dns::RdMatrix HostA  = cla3p::dns::RdMatrix::random(5,5);
    const cla3p::dns::RdVector HostBv = cla3p::dns::RdVector::random(5);
    const cla3p::dns::RdMatrix HostBm = cla3p::dns::RdMatrix::random(5,3);

    /*
     * Transfer to device
     */
    culite::dns::RdMatrix A ;
    culite::dns::RdVector b;
    culite::dns::RdMatrix B;
    HostA >> A;
    HostBv >> b;
    HostBm >> B;

    /*
     * Instantiate LU solver
     */
    culite::LapackLU<culite::dns::RdMatrix> luSolver;

    /*
     * Decompose A into LU product
     */
    luSolver.decompose(A);

    {
        /*
         * Single column (vector) rhs
         * Overwrite x with the solution (A^{-1} * Bv)
         */
        culite::dns::RdVector x = b;
        luSolver.solve(x);
        // Check error norm on host
        cla3p::dns::RdVector HostX;
        x >> HostX;
        std::cout << "Dense Vector rhs::Absolute Error: "
            << (HostBv - HostA * HostX).evaluate().normOne() << std::endl;
    }

    {
        /*
         * Multiple column (matrix) rhs
         * Overwrite X with the solution (A^{-1} * B)
         */
        culite::dns::RdMatrix X = B;
        luSolver.solve(X);
        // Check error norm on host
        cla3p::dns::RdMatrix HostX;
        X >> HostX;
        std::cout << "Dense Matrix rhs::Absolute Error: "
            << (HostBm - HostA * HostX).evaluate().normOne() << std::endl;
    }

    return 0;
}
