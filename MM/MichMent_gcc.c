/*   
**   MichMent_cvode.c
**	 
**   Cvode-Mex implementation of BioNetGen model 'MichMent'.
**
**   Code Adapted from templates provided by Mathworks and Sundials.
**   QUESTIONS about the code generator?  Email justinshogg@gmail.com
**
**   Requires the CVODE libraries:  sundials_cvode and sundials_nvecserial.
**   https://computation.llnl.gov/casc/sundials/main.html
**
**-----------------------------------------------------------------------------
**
**   COMPILE in MATLAB:
**   mex -L<path_to_cvode_libraries> -I<path_to_cvode_includes>  ...
**          -lsundials_nvecserial -lsundials_cvode -lm MichMent_cvode.c
**
**   note1: if cvode is in your library path, you can omit path specifications.
**
**   note2: if linker complains about lib stdc++, try removing "-lstdc++"
**     from the mex configuration file "gccopts.sh".  This should be in the
**     matlab bin folder.
** 
**-----------------------------------------------------------------------------
**
**   EXECUTE in MATLAB:
**   [error_status, species_out, observables_out]
**        = MichMent_cvode( timepoints, species_init, parameters )
**
**   timepoints      : column vector of time points returned by integrator.
**   parameters      : row vector of 4 parameters.
**   species_init    : row vector of 4 initial species populations.
**
**   error_status    : 0 if the integrator exits without error, non-zero otherwise.
**   species_out     : species population trajectories
**                        (columns correspond to states, rows correspond to time).
**   observables_out : observable trajectories
**                        (columns correspond to observables, rows correspond to time).
*/

/* Library headers */
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>             /* prototypes for CVODE  */
#include <nvector/nvector_serial.h>  /* serial N_Vector       */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_spgmr.h>       /* prototype for CVSpgmr */

/* Problem Dimensions */
#define __N_PARAMETERS__   4
#define __N_EXPRESSIONS__  10
#define __N_OBSERVABLES__  4
#define __N_RATELAWS__     3
#define __N_SPECIES__      4

/* core function declarations */
int   check_flag  ( void *flagvalue, char *funcname, int opt );
void  calc_expressions ( N_Vector expressions, double * parameters );
void  calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions );
void  calc_ratelaws    ( N_Vector ratelaws,  N_Vector species, N_Vector expressions, N_Vector observables );
int   calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data );

/* user-defined function declarations */


/* user-defined function definitions  */


/* Calculate expressions */
void
calc_expressions ( N_Vector expressions, double * parameters )
{
    NV_Ith_S(expressions,0) = parameters[0];
    NV_Ith_S(expressions,1) = (602.0*NV_Ith_S(expressions,0));
    NV_Ith_S(expressions,2) = parameters[1];
    NV_Ith_S(expressions,3) = parameters[2];
    NV_Ith_S(expressions,4) = parameters[3];
    NV_Ith_S(expressions,5) = (0.01*NV_Ith_S(expressions,1));
    NV_Ith_S(expressions,6) = (1.0*NV_Ith_S(expressions,1));
    NV_Ith_S(expressions,7) = pow(10.0,NV_Ith_S(expressions,2));
    NV_Ith_S(expressions,8) = pow(10.0,NV_Ith_S(expressions,3));
    NV_Ith_S(expressions,9) = pow(10.0,NV_Ith_S(expressions,4));
   
}

/* Calculate observables */
void
calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions )
{
    NV_Ith_S(observables,0) = NV_Ith_S(species,0);
    NV_Ith_S(observables,1) = NV_Ith_S(species,1);
    NV_Ith_S(observables,2) = NV_Ith_S(species,3);
    NV_Ith_S(observables,3) = NV_Ith_S(species,2);

}

/* Calculate ratelaws */
void
calc_ratelaws ( N_Vector ratelaws, N_Vector species, N_Vector expressions, N_Vector observables )
{  
    NV_Ith_S(ratelaws,0) = NV_Ith_S(expressions,7)*NV_Ith_S(species,0)*NV_Ith_S(species,1);
    NV_Ith_S(ratelaws,1) = NV_Ith_S(expressions,8)*NV_Ith_S(species,2);
    NV_Ith_S(ratelaws,2) = NV_Ith_S(expressions,9)*NV_Ith_S(species,2);

}


/* Calculate species derivatives */
int
calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data )
{
    //int         return_val;
    N_Vector *  temp_data;
    
    N_Vector    expressions;
    N_Vector    observables;
    N_Vector    ratelaws;

    /* cast temp_data */
    temp_data = (N_Vector*)f_data;
     
    /* sget ratelaws Vector */
    expressions = temp_data[0];
    observables = temp_data[1];
    ratelaws    = temp_data[2];
       
    /* calculate observables */
    calc_observables( observables, species, expressions );
    
    /* calculate ratelaws */
    calc_ratelaws( ratelaws, species, expressions, observables );
                        
    /* calculate derivatives */
    NV_Ith_S(Dspecies,0) = -NV_Ith_S(ratelaws,0) +NV_Ith_S(ratelaws,1) +NV_Ith_S(ratelaws,2);
    NV_Ith_S(Dspecies,1) = -NV_Ith_S(ratelaws,0) +NV_Ith_S(ratelaws,1);
    NV_Ith_S(Dspecies,2) = NV_Ith_S(ratelaws,0) -NV_Ith_S(ratelaws,1) -NV_Ith_S(ratelaws,2);
    NV_Ith_S(Dspecies,3) = NV_Ith_S(ratelaws,2);


    return(0);
}

int main()
{
    /* variables */
    double par[4] = {1,-2.77,-1,-2};
    double sp[4] = {6.020000000000e+00,6.020000000000e+02,0,0}; // 6.02 is getting printed as 6.0199999999999996e+00  
    double timepoints[11] = {0.0,2000.0,4000.0,6000.0,8000.0,10000.0,12000.0,14000.0,16000.0,18000.0,20000.0}; 
    /* get number of timepoints */
    size_t n_timepoints = 11;
    double *  parameters = par;
    double *  species_init = sp;
    size_t    i;
    size_t    j;
    /* intermediate data vectors */
    N_Vector  expressions;
    N_Vector  observables;
    N_Vector  ratelaws;
    /* array to hold pointers to data vectors */
    N_Vector  temp_data[3];
    /* CVODE specific variables */
    realtype  reltol;
    realtype  abstol;
    realtype  time;
    N_Vector  species;
    void *    cvode_mem;
    int       flag;   
    /* initialize intermediate data vectors */
    expressions  = NULL;
    expressions = N_VNew_Serial(__N_EXPRESSIONS__);
    if (check_flag((void *)expressions, "N_VNew_Serial", 0))
    {
        return 1;
    }
    observables = NULL;
    observables = N_VNew_Serial(__N_OBSERVABLES__);
    if (check_flag((void *)observables, "N_VNew_Serial", 0))
    {
        N_VDestroy_Serial(expressions);
        return 1;
    }
    ratelaws    = NULL; 
    ratelaws = N_VNew_Serial(__N_RATELAWS__);
    if (check_flag((void *)ratelaws, "N_VNew_Serial", 0))
    {   
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);        
        return 1;
    } 
    /* set up pointers to intermediate data vectors */
    temp_data[0] = expressions;
    temp_data[1] = observables;
    temp_data[2] = ratelaws;
    /* calculate expressions (expressions are constant, so only do this once!) */
    calc_expressions( expressions, parameters );      
    /* SOLVE model equations! */
    species   = NULL;
    cvode_mem = NULL;
    /* Set the scalar relative tolerance */
    reltol = 1e-08;
    abstol = 1e-06;
    /* Create serial vector for Species */
    species = N_VNew_Serial(__N_SPECIES__);
    if (check_flag((void *)species, "N_VNew_Serial", 0))
    {  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);
        return 1;
    }
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   NV_Ith_S(species,i) = species_init[i];   }
    
    /* write initial observables populations into species_out */ 
    calc_observables( observables, species, expressions );  
    /*   Call CVodeCreate to create the solver memory:    
     *   CV_ADAMS or CV_BDF is the linear multistep method
     *   CV_FUNCTIONAL or CV_NEWTON is the nonlinear solver iteration
     *   A pointer to the integrator problem memory is returned and stored in cvode_mem.
     */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                 
        return 1;                        
    }                                  
    /*   Call CVodeInit to initialize the integrator memory:     
     *   cvode_mem is the pointer to the integrator memory returned by CVodeCreate
     *   rhs_func  is the user's right hand side function in y'=f(t,y)
     *   T0        is the initial time
     *   y         is the initial dependent variable vector
     */
    flag = CVodeInit(cvode_mem, calc_species_deriv, timepoints[0], species);
    if (check_flag(&flag, "CVodeInit", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                 
        return 1;                        
    }                                     
    /* Set scalar relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                  
        return 1;                        
    }                                        
    /* pass params to rhs_func */
    flag = CVodeSetUserData(cvode_mem, &temp_data);
    if (check_flag(&flag, "CVodeSetFdata", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                  
        return 1;                        
    }                                     
    /* select linear solver */
    flag = CVDense(cvode_mem, __N_SPECIES__);
    if (check_flag(&flag, "CVDense", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                 
        return 1;                        
    }                                  
    
    flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                 
        return 1;                        
    }                                  
    flag = CVodeSetMaxErrTestFails(cvode_mem, 7);
    if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);              
        return 1;                        
    }                                  
    flag = CVodeSetMaxConvFails(cvode_mem, 10);
    if (check_flag(&flag, "CVodeSetMaxConvFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                
        return 1;                        
    }                                  
    flag = CVodeSetMaxStep(cvode_mem, 0.0);
    if (check_flag(&flag, "CVodeSetMaxStep", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);                 
        return 1;                        
    }                                  
    printf("%.16e\t",timepoints[0]);
    for(j=0;j<__N_SPECIES__;j++)
    {
        printf("%.16e\t",NV_Ith_S(species,j));
    }
    printf("\n");
    /* integrate to each timepoint */
    for ( i=1;  i < n_timepoints;  i++ )
    {
        flag = CVode(cvode_mem, timepoints[i], species, &time, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1))
        {
            N_VDestroy_Serial(expressions);
            N_VDestroy_Serial(observables);           
            N_VDestroy_Serial(ratelaws);
            N_VDestroy_Serial(species);
            CVodeFree(&cvode_mem);
            return 1;
        } 
        printf("%.16e\t",timepoints[i]);
        for(j=0;j<__N_SPECIES__;j++)
        {
            printf("%.16e\t",NV_Ith_S(species,j));
        }
        printf("\n");
    }
    /* Free vectors */
    N_VDestroy_Serial(expressions);
    N_VDestroy_Serial(observables);  
    N_VDestroy_Serial(ratelaws);        
    N_VDestroy_Serial(species);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return flag;
}
/*  Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        printf( "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname );    
        return(1);
    }
    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0)
        {
            printf( "\nSUNDIALS_ERROR: %s() failed with flag = %d\n", funcname, *errflag );
            return(1);
        }
    }
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        printf( "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n", funcname );
        return(1);
    }
    return(0);
}
