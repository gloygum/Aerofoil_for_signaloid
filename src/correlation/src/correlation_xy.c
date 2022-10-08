/*
    Example to investigate correlation between two random variables, one of which
    is linearly dependent on the other, using  Signaloid unvertain lib.
    Clive Emary 08-10-22
*/

#include <math.h>
#include <stdio.h>
#include <uncertain.h>



// input loader
static void loadInputs(double *  , double *  );

// Pearson Correlation coefficient
static double correlation(double , double );


int
main(int argc, char *	argv[])
{
   	// x,z are indepdent r.v.
	double	x,z;

	loadInputs(&x, &z);

   	// y in a linear combination of x and z
	double y = x+z;
	
	printf("x and z and independent random vairables. y =(x+z)/2\n\n");
	
	printf("1st method \n");

   	// examine distribution of y
	printf("Distribution of y:\n");
	printf("y = %.2E \n",y);
	libUncertainDoublePrint(y);

	// CE comment
	printf("This distribution is correct.\n");
	printf("\nCorrelation: rho(x,y) = %.2f\n",correlation(x,y));
	printf("This correlation is not correct. It should be 1/Sqrt[2]\n");

	// Alternative method without input loader
	// in case this makes a difference
	x  = libUncertainDoubleGaussDist(1, 1 );
	z  = libUncertainDoubleGaussDist(1, 1 );
	y = x+z;

	// CE comment
	printf("\n\n2nd method \n");
	printf("y = %.2E \n",y);
	libUncertainDoublePrint(y);
	printf("Correlation rho(x,y) = %.2f\n",correlation(x,y));
	
	printf("Also not correct \n");




	// evaluate the pieces by hand
	printf("If I decompose the correlation function by hand and evaluate the bits individually, I get\n");
	// numerator = <xy> - <x><y>
	double numerator = 0.5*libUncertainDoubleNthMoment(x*x,1) + 0.5*libUncertainDoubleNthMoment(x*z,1)
		- libUncertainDoubleNthMoment(x,1)*0.5*(libUncertainDoubleNthMoment(x,1)+libUncertainDoubleNthMoment(z,1));

	// denominator = \sqrt(<x^2>_c <y^2>_c )
	double denominator = sqrt(
				libUncertainDoubleNthMoment(x,2) * 0.25 * (libUncertainDoubleNthMoment(x,2) + libUncertainDoubleNthMoment(z,2))
			    );
	//CE comment
	printf("   rho(x,y) = %.02f\n",numerator/denominator);
	printf("This gives the correct answer when run on a core with autocorrelation tracking....");

	return 0;
}



static void
loadInputs(double *  x, double *  z)
{
	*x  = libUncertainDoubleGaussDist(1, 1 );
	*z  = libUncertainDoubleGaussDist(1, 1 );

}

/* Correlation coefficient
    Inputs:
            X , Y: double random variables with attached signaloid distribution
    Outputs:
            Pearson correlation coefficient

            rho(X,Y) = ( <XY>-<X><Y> ) / (sigma_X sigma_Y)

            with sigma_{X,Y} the standard deviations of the two distributions
*/
static double correlation(double X, double Y){
	double rho = ( libUncertainDoubleNthMoment(X*Y,1) -  libUncertainDoubleNthMoment(X,1)*  libUncertainDoubleNthMoment(Y,1))
		 / sqrt(libUncertainDoubleNthMoment(X,2) *  libUncertainDoubleNthMoment(Y,2));

	return rho;
}
