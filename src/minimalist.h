
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

class CurvatureSearchScope{
	unsigned int size;
	double* lastderiv;
	double lastvalue;
	double expdiff;
	double* curvature;
	double log_alpha; // temp... to remove
    unsigned int shrinkcount;
public:

	CurvatureSearchScope(): lastderiv(NULL), curvature(NULL){}
	CurvatureSearchScope(const CurvatureSearchScope&);
	~CurvatureSearchScope() {delete[](lastderiv);}
	CurvatureSearchScope& operator=(const CurvatureSearchScope&);

	void init(double initalalpha, unsigned int in_size);

	double updateAscent(const double &value, double * guess,  double const * const deriv); // returns a safe bound on the parameter changes;
	double updateDescent(const double &value, double * guess , double const * const deriv); // returns a safe bound on the parameter changes;

	double checkDerivative(const double &value, double * guess , double const * const deriv); // test if input derivative are correct

	void wrFinalGuess(double* guess)const;
    double getLastValue()const{return lastvalue;}
    void show(FILE* f = stdout, int level=0) const{fprintf(f,"Curvature scope in %i dimentions, last value %e\t expected difference %e, surrent log-alpha %e\n",size,lastvalue,expdiff,log_alpha);}
};
