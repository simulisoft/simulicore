/**
 * @example ex01b_dense_vector_fill.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
    cla3p::dns::RdVector x(3);
    cla3p::dns::RdVector y(3);
    cla3p::dns::RdVector z(3);

    /*
     * Set all values in x to 1.5 (using member function).
     * Set all values in y to -2.7 (using operator=).
     */
    x.fill(1.5);
    y = -2.7;

    /*
     * Insert values in z manually.
     */
    z(0) = 3.;
    z(1) = 2.;
    z(2) = z(0) + z(1);

    std::cout << "x:\n" << x << "\n";
    std::cout << "y:\n" << y << "\n";
    std::cout << "z:\n" << z;

    return 0;
}
