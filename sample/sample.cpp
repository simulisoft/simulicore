#include <iostream>
#include <omp.h>

#include <cla3p/dense.hpp>

#if !defined(__APPLE__)
#include <culite/dense.hpp>
#endif

/*-------------------------------------------------*/
static void test_omp()
{
    #pragma omp parallel
	{
		#pragma omp critical
		std::cout << "Hello from thread " << omp_get_thread_num() << " / " << omp_get_max_threads() << std::endl;
	}
}
/*-------------------------------------------------*/
static void test_cla3p()
{
    std::cout << "CLA3P Integer size is: " << sizeof(cla3p::int_t) << std::endl;

    cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(5,3);
    std::cout << A.info("A") << A << std::endl;

    cla3p::dns::RdMatrix At = A.transpose();
    std::cout << At.info("A transpose") << At << std::endl;
}
/*-------------------------------------------------*/
static void test_culite()
{
#if !defined(__APPLE__)
    std::cout << "cuLite Integer size is: " << sizeof(culite::int_t) << std::endl;

    cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(5,3);
    std::cout << A.info("host A") << A << std::endl;

    culite::dns::RdMatrix devA;
    A >> devA;
    std::cout << devA.info("device A") << devA << std::endl;

    culite::dns::RdMatrix devAt = devA.transpose();
    std::cout << devAt.info("device A transpose") << devAt << std::endl;
#else
    std::cout << "culite is not supported on Apple platforms due to lack of CUDA support." << std::endl;
#endif
}
/*-------------------------------------------------*/
static void welcome()
{
    std::cout << "Welcome to SimuliCore!" << std::endl;
}
/*-------------------------------------------------*/
int main()
{
    welcome();
    test_omp();
    test_cla3p();
    test_culite();

	return 0;
}
/*-------------------------------------------------*/
