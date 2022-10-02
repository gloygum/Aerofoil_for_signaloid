/*
 * Aerofoil lift with Signaloid uncertainties
 * Clive Emary 02-10-2022
 */

/*
 * Calculates Mach number and the lift of an aerofoil from Pitot tube pressure measurements.
 * The equations used here are valid for subsonic, compressible flow with
 * a thin aerofoil with low angle of attack.
 * See Anderson, J. D. Jr., Fundamentals of Aerodynamics, Sixth Edition McGraw-Hill (2017).
 *
 * Inputs
 *      pressure_fraction : excess pressure fraction
 *      gamma :             specific heat ratio of air
 *      c :                 chord length of aerofoil
 *      alpha :             angle of attack
 *      P_static :          static pressure
 *      P_pitot :           presure in Pitot tube
 *
 * Outputs
 *      mach :              mach number
 *      lift :              lift per unit span
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <uncertain.h>


// PARAMETERS

// fix % to be used for the evaluation of confidence intervals
const double level = 95.;

// verbose = 1 gives full output to stdout; verbose = 0, just the minimum
const int verbose = 0;

// Output file
const char * fname = "liftdata/lift.csv";



// FUNCTION DECLARATIONS

// INPUT LOADER
static void loadInputs(double, double * , double * , double * , double * , double * );

// CONFIDENCE INTERVAL
static double confidence_interval(double , double );

// CORRELATION COEFFICIENT
static double correlation(double , double );



// MAIN
int
main(int argc, char *	argv[])
{
	double	gamma, c, alpha, P_static, P_pitot;
    double  mach, lift;
    double mean_M, mean_L, Delta_M, Delta_L;

    double frac;
    double p_frac_example = 0.1;

    // open data file & exit if this failed
    FILE * fp = fopen(fname, "w");
    if (fp == NULL)
    {
        printf("Could not open file. Check you've mounted a data source");
        return 0;
    }


    //  Loop over values of pressure_fraction
    printf("Looping over pressure fraction...\n\n");

    for (frac = 0.01; frac < 0.9; frac=frac+0.05)
    {

	    loadInputs(frac, &gamma, &c, &alpha, &P_static, &P_pitot);

        // should ASCII the function into here
        mach = sqrt(2.0/(gamma - 1.0)
               *  (pow(fmax(P_pitot/P_static,1.0), 1.0 - (1.0/gamma)) - 1.0));

        // should ASCII the function into here
        lift = c * M_PI * alpha * gamma * P_static * pow(mach, 2.0)
               / sqrt(1.0- pow(mach,2.0));

        // summary statistics
        mean_M = libUncertainDoubleNthMoment(mach,1);
        Delta_M = confidence_interval(mach,level);

        mean_L = libUncertainDoubleNthMoment(lift,1);
        Delta_L = confidence_interval(lift,level);


        // report if verbose  = 1
        if (verbose == 1)
        {
            printf("pressure fraction %.3f\n", frac);

            printf("Mach number mean & %.0fpct interval: %.2f, [%.2f. %.2f]\n",
                   level, mean_M,  mean_M - Delta_M, mean_M + Delta_M );

            printf("Lift/unit spanmean & %.0fpct interval: %.2f, [%.2f. %.2f]\n\n",
                   level,mean_L,mean_L - Delta_L, mean_L + Delta_L );
        }

        // Write data to file
        fprintf(fp,"%f, %f, %f, %f, %f, %f, %f\n",
                frac,
                mean_M, mean_M - Delta_M, mean_M + Delta_M,
                mean_L, mean_L - Delta_L, mean_L + Delta_L);

    }

    // close file to finish loop
    fclose(fp);
    printf("...Loop complete.\n\n");
    printf(
        "Mach & Lift data writen to '%s' as CSV.\n"
        "You can tidy it up with:\n"
        "a sed command that somehow wrecks the output when verbose = 1 "
        "but is okay when verbose = 0 ...\n"
        // uncomment following line to see what happens ...
        // "sed -E 's/Ux[[:alnum:]]*//g' lift.csv | tee lift2.csv\n"
        ,
        fname
    );


    // single example for discussion
    printf("\n------------------------------\n");
    printf("Single example for discussion\n");
    printf("------------------------------\n");

    frac = p_frac_example;
    printf("pressure fraction %.3f\n", frac);

	loadInputs(frac, &gamma, &c, &alpha, &P_static, &P_pitot);

    // should ASCII the function into here
    mach = sqrt(2.0/(gamma - 1.0)
           *  (pow(fmax(P_pitot/P_static,1.0), 1.0 - (1.0/gamma)) - 1.0));

    // should ASCII the function into here
    lift = c * M_PI * alpha * gamma * P_static * pow(mach, 2.0)
           / sqrt(1.0- pow(mach,2.0));



    // summary statistics
    mean_M  = libUncertainDoubleNthMoment(mach,1);
    Delta_M = confidence_interval(mach,level);

    mean_L  = libUncertainDoubleNthMoment(lift,1);
    Delta_L = confidence_interval(lift,level);


    // report distributions
    printf("\nMach number %.02f\n", mach);
    printf("Lift/span %.02E\n\n", lift);

    //report derived quantities
    printf("Mach number mean & %.0fpct interval: %.2f, [%.2f. %.2f]\n",
           level, mean_M,  mean_M - Delta_M, mean_M + Delta_M);
    printf("Lift/unit span mean & %.0fpct interval: %.2f, [%.2f. %.2f]\n\n",
           level,mean_L,mean_L - Delta_L, mean_L + Delta_L);


    //report pressure-pressure correlation function
    printf("pressure-pressure correlation, rho = %.2f.\n",
            correlation(P_pitot, P_static)
    );
    printf("Should be approx 0.7 for s=1/2 ...");



	return 0;
}



// FUNCTION DEFINTIONS


// INPUT LOADER
static void loadInputs(double pressure_fraction,
                       double *  gamma,
                       double *  c,
                       double *  alpha,
                       double *  P_static,
                       double *  P_pitot)
{
    //  input model parameters
    double mean_p_static = 1.01325E5;
    double sigma_p = 0.01 * mean_p_static;
    double sigma_det = 0.0015 *  mean_p_static;

    double mean_alpha = 4.;
    double sigma_alpha = 0.08;

    double mean_c = 0.305;
    double tol_c = 4E-4;

    double s = 0.5;
    //  include gamma here like this; temperature & in particular humdity variations
    //  could result in a distribution for gamma
    *gamma = 1.400;

    //  maximum pressure in Pitot tube at Mach = 1
    double max_p_pitot = mean_p_static
                         * pow( (*gamma-1.0)/2.0 * 1.0 + 1.0, *gamma/(*gamma - 1.0));

    //  mean pitot pressure for further calculations given as presure_fraction
    //  of max value above statics pressure
    double mean_p_pitot = mean_p_static + pressure_fraction
                              * (max_p_pitot -  mean_p_static);

    //  distributions
    double XP = libUncertainDoubleGaussDist( mean_p_static, sigma_p);

    *P_static = XP + libUncertainDoubleGaussDist(0.0, sigma_det);

    *P_pitot  = s * (XP - mean_p_static + mean_p_pitot)
                + (1-s) * libUncertainDoubleGaussDist(mean_p_pitot,sigma_p)
                + libUncertainDoubleGaussDist(0.0, sigma_det);

                *alpha    = libUncertainDoubleGaussDist( mean_alpha, sigma_alpha );

    *c        = libUncertainDoubleUniformDist( mean_c-tol_c, mean_c+tol_c );

}


/*  CONFIDENCE INTERVAL
 *
 *  Inputs:
 *      dist : double variable with attached signaloid distribution
 *      level : size of confidence interval, specified as percent
 *  Output:
 *      Delta : symmetric width of confidence interval from mean
 *              calculated from P( -Delta < X- mean(x) < Delta ) = level/100
 *                where P(X) is the distribution in question.
*/
static double confidence_interval(double dist, double level)
{
    double mean =  libUncertainDoubleNthMoment(dist,1);
    double step = (libUncertainDoubleSupportMax(dist)
                     - libUncertainDoubleSupportMin(dist))/10000.;
    double Delta = 0.0;
    double cum_prob = 0.0;

    while(cum_prob < level/100.)
    {
        Delta = Delta + step;
        cum_prob = libUncertainDoubleProbabilityGT(dist,  mean - Delta)
                       - libUncertainDoubleProbabilityGT(dist,  mean + Delta);
    }

    return Delta;
}


/* CORRELATION COEFFICIENT
 *
 * Inputs:
 *          X, Y :  double random variables with attached signaloid distribution
 * Outputs:
 *          rho :   Pearson correlation coefficient
 *                  rho(X,Y) = ( <XY>-<X><Y> ) / (sigma_X sigma_Y)
 *                  with sigma_{X,Y} the standard deviations of the two distributions
*/
static double correlation(double X, double Y)
{
    double rho = (libUncertainDoubleNthMoment(X*Y,1)
                       - libUncertainDoubleNthMoment(X,1)*  libUncertainDoubleNthMoment(Y,1))
                 / sqrt(libUncertainDoubleNthMoment(X,2) *  libUncertainDoubleNthMoment(Y,2));

    return rho;
}
