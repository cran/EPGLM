
#include <boost/random/exponential_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/random/student_t_distribution.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/random/binomial_distribution.hpp>
#define _mINF -100000000000
#define tZERO 0.000000000001
#define _THRES 0.0000001
#define MAXIT 200
#define NSTEP 3
#define MAXBIT 30
#define SCRAM 3
#define INIT 1 



//standard template library
#include <string>
#include <algorithm>
#include <vector>
#include <list>
#include <math.h>


using namespace boost::random;
using namespace boost::math;
using namespace arma;
using namespace std;
using namespace boost::math;

#include "EP.h"
#include "EPclogit.h"
