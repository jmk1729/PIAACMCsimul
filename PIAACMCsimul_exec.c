/**
 * @file    PIAAACMCsimul_init.c
 * @brief   PIAA-type coronagraph design, execute
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    21 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/stat.h>
#include <time.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern DATA data;   

extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;








/** @brief Make Lyot stops using off-axis light minimums
 * 
 *  Finds minumum flux level in 3D intensity data cube
 * 
 * @param[in]  IDincohc_name             image: Input 3D intensity 
 * @param[out] test_oals_val.fits        FITS : Minimum 2D image
 * @param[out] test_oals_index.fits      FITS : z-index for each pixel of the 2D output minimum
 * 
*/

long PIAACMCsimul_optimizeLyotStop_OAmin(const char *IDincohc_name)
{
    long IDincohc;
    
    long IDindex;
    long IDminflux;
    long ii, jj, kk;
    long xsize = 0;
    long ysize = 0;
    long xysize = 0;
    double minv;
    long minindex;
    double tmpv;
    long NBz;
    

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif

    
    IDincohc = image_ID(IDincohc_name);

    xsize = data.image[IDincohc].md[0].size[0];
    ysize = data.image[IDincohc].md[0].size[1];
    xysize = xsize*ysize;
    NBz = data.image[IDincohc].md[0].size[2];
    
    IDindex = create_2Dimage_ID("oals_index", xsize, ysize);
    IDminflux = create_2Dimage_ID("oals_val", xsize, ysize);
    
    for(ii=0;ii<xsize;ii++)
        for(jj=0;jj<ysize;jj++)
            {
                minv = data.image[IDincohc].array.F[jj*xsize+ii];
                minindex = 0;

                for(kk=1;kk<NBz;kk++)
                {
                    tmpv = data.image[IDincohc].array.F[kk*xysize+jj*xsize+ii];
                    if(tmpv<minv)
                        {
                            minv = tmpv;
                            minindex = kk;
                        }
                }
                data.image[IDminflux].array.F[jj*xsize+ii] = minv;
                data.image[IDindex].array.F[jj*xsize+ii] = minindex;
            }    
    
    save_fits("oals_index", "!test_oals_index.fits");
    save_fits("oals_val", "!test_oals_val.fits");
    
    return 0;
}


























		/** 
		 * ---
		 * 
		 * ## Mode 3: Calibrate, no focal plane mask (not currently used)
		 * 
		 * Compute PSF and contrast with no focal plane mask with the current design.
        * 
        * Provides the denominator for the contrast estimate
        * 
        * Saved by PIAACMCsimul_computePSF as fits file "psfi0"
        
        */
double PIAACMCsimul_exec_mode03()
{   
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;    
	double val;
    double paramref[10000];
    
	        printf("=================================== mode 003 ===================================\n");

       // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }


        // init as in mode 0
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        paramref[0] = piaacmc[0].fpmaskamptransm;

        piaacmc[0].fpmaskamptransm = -1.0;  // Remove focal plane mask
        piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
        // compute the PSF for an on-axis source, all optical elements
        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

        // restore original configuration
        piaacmc[0].fpmaskamptransm = paramref[0];
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
        PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);
        piaacmcsimul_var.FORCE_CREATE_fpmza = 0;

return val;
}






    /** 
     * ---
     * 
     * ## Mode 4: Optimize PIAA optics shapes, cosine modes only (not currently used, replaced by mode 40. skipping)
     * 
     */ 
int PIAACMCsimul_exec_mode04()
{   
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;      
    long NBiter = 1000;    
    long number_param;
    long kmax;
    long k;    
    
    int paramtype[10000]; // FLOAT or DOUBLE
    float *paramvalf[10000]; // array of pointers, float
    double paramrefval[10000];
 
    double paramdelta[10000];
    double parammaxstep[10000]; // maximum single iteration step
    double parammin[10000]; // minimum value
    double parammax[10000]; // maximum value
 
            printf("=================================== mode 004 ===================================\n");
       // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        piaacmcsimul_var.LINOPT = 1; // perform linear optimization
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value.f+0.01;
        else
            NBiter = 1000;

        kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmax = (long) data.variable[IDv].value.f+0.01;

        if(kmax>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];

        number_param = 0;
        for(k=0; k<kmax; k++)
        {
            paramtype[number_param] = _DATATYPE_FLOAT;
            paramvalf[number_param] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            paramdelta[number_param] = 1.0e-9;
            parammaxstep[number_param] = 1.0e-8;
            parammin[number_param] = -1.0e-5;
            parammax[number_param] = 1.0e-5;
            number_param++;
        }

        for(k=0; k<kmax; k++)
        {
            paramtype[number_param] = _DATATYPE_FLOAT;
            paramvalf[number_param] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            paramdelta[number_param] = 1.0e-9;
            parammaxstep[number_param] = 1.0e-8;
            parammin[number_param] = -1.0e-5;
            parammax[number_param] = 1.0e-5;
            number_param++;
        }
        piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 1;
        piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 1;


return 0;
}









   /** 
     * ---
     * 
     * ## Mode 5: Optimize Lyot stops shapes and positions
     * 
     * 
     * 
     */ 
int PIAACMCsimul_exec_mode05()
{
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;   
	long elem, elem0;
    long ls, NBpropstep, k;
    double oaoffset; // off-axis offset for Lyot stop optimization	
    double zmin, zmax;	
    char fnamea[500];
    char fnamep[500];	
    double lstransm;
	long NBincpt, k1;
	long ii;
    long xsize = 0;
    long ysize = 0;
    long NBkr, kr;
    long ID1, IDa, IDc;
    char command[1000];   
    char fname[1000];  
    long cnt;    
    char fptestname[1000];
    FILE *fptest;
    int r;
     	
        printf("=================================== mode 005 ===================================\n");
       // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        
        
        /// ### Initialize as in mode 0
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm



        /// ### Load CLI variables as appropriate
             
		/// - <- **PIAACMC_nbpropstep** : # of propagation steps along the beam
        NBpropstep = 150;
        if((IDv=variable_ID("PIAACMC_nbpropstep"))!=-1)
            NBpropstep = (long) data.variable[IDv].value.f+0.01;
              
		/// - <- **PIAACMC_lstransm** : desired Lyot stop transmission		
        lstransm = 0.85;
        if((IDv=variable_ID("PIAACMC_lstransm"))!=-1)
            lstransm = (double) data.variable[IDv].value.f;
        printf("lstransm  = %f\n", lstransm);





        /// ### Identify post focal plane pupil plane (first pupil after focal plane mask)

        /// Provides reference complex amplitude plane for downstream analysis
        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
        
        printf("=========== %ld elements ======================================================\n", optsyst[0].NBelem);
        // find the ID of the "post focal plane mask pupil" element
        for(elem=0; elem<optsyst[0].NBelem; elem++)
        {
            if( strcmp("post focal plane mask pupil", optsyst[0].name[elem]) == 0 )
            {
                elem0 = elem;
                printf("post focal plane mask pupil = %ld\n", elem);
            }
            else
                printf("elem %ld : %s\n", elem, optsyst[0].name[elem]);
        }
        optsyst[0].keepMem[elem0] = 1; // keep it for future use
        
        
        
        
        
		/// ### Compute incoherent 3D illumination near pupil plane
		
		/// Multiple off-axis sources are propagated and the corresponding intensities added
		/// - -> **OAincohc**  Output incoherent image (3D)
		
        oaoffset = 20.0; // off axis amplitude
        // compute the reference on-axis PSF
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

        // filenames of the complex amplitude and phase in the post FPM pupil plane indexed by elem0
        sprintf(fnamea, "WFamp0_%03ld", elem0);
        sprintf(fnamep, "WFpha0_%03ld", elem0);

        printf("elem0 = %ld\n", elem0);

        // args for the PIAACMCsimul_CA2propCubeInt function
       
        // set range of propagation
        zmin = piaacmc[0].LyotZmin;
        zmax = piaacmc[0].LyotZmax;
            
        // args that determine "extended" off-axis source
        // number of off-axis sources on each circle
        NBincpt = 15;
        // number of circle radii
        NBkr = 5;

        // propagate complex amplitude in a range from zmin to zmax, where 0 is elem0
        // computes the diffracted light from the on-axis source
        ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
        // complex amplitude at elem0, only used to determine image size
        IDa = image_ID(fnamea);

        xsize = data.image[IDa].md[0].size[0];
        ysize = data.image[IDa].md[0].size[1];

        /// OAincohc is the summed light "all" from off-axis sources in the pupil,
        /// including the on-axis source(!),
        /// giving intensity contribution of all off-axis sources
        /// in order to preserve the intensity of the off-axis in the design.
        /// load OAincohc if exist, maybe we've been here before
        sprintf(fname, "%s/OAincohc.fits", piaacmcsimul_var.piaacmcconfdir);
        IDc = load_fits(fname, "OAincohc", 1);


        if(IDc==-1) // OAincohc does not exist so we have to make it
        {
            // create image to receive sum
            IDc = create_3Dimage_ID("OAincohc", xsize, ysize, NBpropstep);

            cnt = 0; // initialize counter so later we can normalize by number of sources
            // loop over radii
            for(kr=0; kr<NBkr; kr++)
            {
                // loop over points at current radius
                for(k1=0; k1<NBincpt; k1++)
                {
                    // compute PSF for a point at this angle with scaled offset
                    // PIAACMCsimul_computePSF changes fnamea and fnamep (in call to OptSystProp_run)!
                    PIAACMCsimul_computePSF(oaoffset*(1.0+kr)/NBkr*cos(2.0*M_PI*k1/NBincpt), oaoffset*(1.0+kr)/NBkr*sin(2.0*M_PI*k1/NBincpt), 0, optsyst[0].NBelem, 0, 0, 0, 0);
                    // propagate that elem0 from zmin to zmax with new PSF
                    ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
                    for(ii=0; ii<xsize*ysize; ii++) // ii is indexing x-y plane
                        for(k=0; k<NBpropstep; k++) // k is indexing z-direction. adding to IDc which looks right
                            data.image[IDc].array.F[k*xsize*ysize+ii] += data.image[ID1].array.F[k*xsize*ysize+ii];
                    delete_image_ID("iproptmp");
                    cnt ++;
                }
            }
            // scale by the number of sources to give average
            for(ii=0; ii<xsize*ysize; ii++)
                for(k=0; k<NBpropstep; k++)
                    data.image[IDc].array.F[k*xsize*ysize+ii] /= cnt;


            sprintf(fname, "!%s/OAincohc.fits", piaacmcsimul_var.piaacmcconfdir);
            // save the final result
            save_fits("OAincohc", fname);
        }
        
        
        
        /// ### Compute on-axis PSF 3D intensity to define light to reject
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        // propagate it into the optical system, with result in image named "iprop00"
        
        /// call function PIAACMCsimul_CA2propCubeInt() to compute 3D intensity cube
        /// @param[out] iprop00 3D intensity image
        PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iprop00");
        //  save_fits("iprop00", "!test_iprop00.fits");

        /// ### Compute image that has the min along z of OAincohc at each x,y
		/// Function PIAACMCsimul_optimizeLyotStop_OAmin() computes minimal intensity image.
		/// 
        PIAACMCsimul_optimizeLyotStop_OAmin("OAincohc");
        
        /// ### Optimize Lyot stops
        
        /// The actual Lyot stop shape and location optimization is done by PIAACMCsimul_optimizeLyotStop(), producing optimal Lyot stops in optLM*.fits
        /// and position relative to elem0 in piaacmc[0].LyotStop_zpos
        PIAACMCsimul_optimizeLyotStop(fnamea, fnamep, "OAincohc", zmin, zmax, lstransm, NBpropstep, piaacmc[0].NBLyotStop);

        sprintf(fptestname, "conj_test.txt");
        fptest = fopen(fptestname, "w");
        fprintf(fptest, "# Optimal Lyot stop conjugations\n");
        fprintf(fptest, "# \n");
        fprintf(fptest, "# DO NOT EDIT THIS FILE\n");
        fprintf(fptest, "# Written by %s in %s\n", __FUNCTION__, __FILE__);
		fprintf(fptest, "# \n");
        fprintf(fptest, "# Lyot stop index   zmin   zmax    LyotStop_zpos    elemZpos[elem0]\n");
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            fprintf(fptest, "%5ld  %f  %f     %f  %f\n", ls, zmin, zmax, piaacmc[0].LyotStop_zpos[ls], optsyst[0].elemZpos[elem0]);
        fclose(fptest);

        // convert Lyot stop position from relative to elem0 to absolute
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] += optsyst[0].elemZpos[elem0];

        // and we're done!  save.
        PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);
        // copy to the final Lyot stop file for this mode
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
        {
            sprintf(command, "cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits", piaacmcsimul_var.piaacmcconfdir, ls, piaacmcsimul_var.piaacmcconfdir, ls);
            r = system(command);
        }

return 0;
}











    /** 
     * ---
     * 
     * ## Mode 11: Setup multizone ring mask and Compute polychromatic response to zones, store result in FPMresp
      
       here we compute how the light propagates from each individual mask zone to the focal plane
       (where each mask zone is completely tranparent)
       * 
       */ 


int PIAACMCsimul_exec_mode11()
{
	long IDv;
    int tmpnblambda;	
    long tmpNBrings;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;  	
	long index;
    char fname[1000];     
    char fnamet[1000]; 
    char fnametmp[1000];  
    char fname1[500];
    char fname2[500];   
    long ID;
	long elem, elem0;

    char command[1000]; 

    // spreading computations over multiple processes for resp matrix
    long mzoffset = 0;
    long mzstep = 1;
    long thr;
    long tmpl1;
    
    long IDcomb;
    double val;
    
	long i, ii, k;
	FILE *fpt;

	int r;
	
	char fnamecomb[1000];
	long ID1;
	
	long IDfpmresp;
    long mz;


	        printf("=================================== mode 011 ===================================\n");

        // get cli variables

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        


        if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
            tmpnblambda = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_NBrings"))!=-1)
            tmpNBrings = data.variable[IDv].value.f;

        // initialize
        PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);

        printf("piaacmcconfdir     : %s\n", piaacmcsimul_var.piaacmcconfdir);
        fflush(stdout);
        printf("SCORINGMASKTYPE    : %d\n", piaacmcsimul_var.SCORINGMASKTYPE);
        fflush(stdout);
        printf("PIAACMC_FPMsectors : %d\n", piaacmcsimul_var.PIAACMC_FPMsectors);
        fflush(stdout);
        printf("lamda              : %ld nm\n", (long) (1.0e9*piaacmc[0].lambda + 0.1));
        fflush(stdout);
        printf("lamdaB             : %ld \n", (long) (1.0*piaacmc[0].lambdaB + 0.1));
        fflush(stdout);
        printf("piaacmc[0].NBrings : %ld\n", piaacmc[0].NBrings);
        fflush(stdout);
        printf("mask rad           : %ld\n", (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1));
        fflush(stdout);
        printf("computePSF_ResolvedTarget : %d\n", piaacmcsimul_var.computePSF_ResolvedTarget);
        fflush(stdout);
        printf("computePSF_ResolvedTarget_mode : %d\n", piaacmcsimul_var.computePSF_ResolvedTarget_mode);
        fflush(stdout);
        printf("piaacmc[0].fpmmaterial_name : %s\n", piaacmc[0].fpmmaterial_name);
        fflush(stdout);
        printf("piaacmc[0].nblambda         : %d\n", piaacmc[0].nblambda);
        fflush(stdout);

        // set output filename of the combined focal plane mask response file
		sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);
  
        printf("fname = %s\n", fname);
        fflush(stdout);


        // get the combined focal plane mask response
        ID = load_fits(fname, "FPMresp", 1);
        // if it did not exist, create it
        if(ID==-1)
        {

            // get the number of tmux threads from cli
            piaacmcsimul_var.PIAACMC_FPMresp_mp = 1; // 1: all computations on a single thread
            if((IDv=variable_ID("PIAACMC_FPMresp_mp"))!=-1) // multi threaded
                piaacmcsimul_var.PIAACMC_FPMresp_mp = (long) data.variable[IDv].value.f+0.01;
            printf("PIAACMC_FPMresp_mp = %ld\n", piaacmcsimul_var.PIAACMC_FPMresp_mp);

            printf("------------------------------------- STEP02\n");
            printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            fflush(stdout);

            // get our tmux thread number in [0 PIAACMC_FPMresp_mp]
            // where the master thread has PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp
            piaacmcsimul_var.PIAACMC_FPMresp_thread = 0;
            if((IDv=variable_ID("PIAACMC_FPMresp_thread"))!=-1) // multi threaded
                piaacmcsimul_var.PIAACMC_FPMresp_thread = (long) data.variable[IDv].value.f+0.01;
            printf("PIAACMC_FPMresp_thread = %ld\n", piaacmcsimul_var.PIAACMC_FPMresp_thread);


            
            index = 0;
            if((piaacmcsimul_var.PIAACMC_FPMresp_mp==1)||(piaacmcsimul_var.PIAACMC_FPMresp_thread > piaacmcsimul_var.PIAACMC_FPMresp_mp-1))  // main or combine process
                                            // why not test PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp?
            {
                // we're the parent set up the FPM zone map
                piaacmcsimul_var.FORCE_CREATE_fpmzmap = 1;
                piaacmcsimul_var.FORCE_CREATE_fpmzt = 1;
                piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
                PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
            }
            else
            {
                printf("NO initOK file created\n");
                // we're a child tmux thread, do not set up the FPM zone map, get it from parent via file
                piaacmcsimul_var.FORCE_CREATE_fpmzmap = 0;
                piaacmcsimul_var.FORCE_CREATE_fpmzt = 0;
                piaacmcsimul_var.FORCE_CREATE_fpmza = 0;
                PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
            }



            // if we're the parent load
            if((piaacmcsimul_var.PIAACMC_FPMresp_mp==1)||(piaacmcsimul_var.PIAACMC_FPMresp_thread > piaacmcsimul_var.PIAACMC_FPMresp_mp-1))
            {
                sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);

                sprintf(fname1, "!%s.tmp", fname);
                sprintf(fname2, "!%s", fname);
                mzoffset = 0;
                mzstep = 1;

                ID = load_fits(fname, "FPMresp", 1);  // this will always fail in the current state (see line 6606) ************************
                IDcomb = ID;
            }
            else // we're a child tmux thread.
            {
                // combined FPMresp file
               sprintf(fnamecomb, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);


                // partial FPMresp file
               sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d_mp%02ld_thread%02ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda, piaacmcsimul_var.PIAACMC_FPMresp_mp, piaacmcsimul_var.PIAACMC_FPMresp_thread);

                // stash the filename of the partial file for later
                sprintf(fname1, "!%s.tmp", fname);
                sprintf(fname2, "%s", fname);
                // set region of the partial file that this child computes
                mzoffset = piaacmcsimul_var.PIAACMC_FPMresp_thread;
                mzstep = piaacmcsimul_var.PIAACMC_FPMresp_mp;

                ID = load_fits(fname, "FPMresp", 1); // may exist from a previous execution with restart
                IDcomb = load_fits(fnamecomb, "FPMresp", 1); // will always fail
            }
            // at this point IDcomb==-1, and in the parent ID==-1 always, and in the child ID==-1 if this is not a restart
            // actually create the FPMresp file either as a part by a child or combined by the parent
            if((IDcomb==-1)&&(ID==-1))  // this will always fire for the parent thread,
                                            // and will always fire for children in a fresh run
            {
                //                printf("--------------------------------------------------------STEP 0005 File \"%s\" does not exist: creating\n", fname);
                //   printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //  sleep(3);

                // usual initialzation
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                //piaacmc[0].fpmaskamptransm = 1.0;
                // set the physical size of the FPM as mean(lambda/D)*mask radius in units of lambda/D
                piaacmc[0].fpmRad = 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio * piaacmcsimul_var.PIAACMC_MASKRADLD; // piaacmcsimul_var.PIAACMC_MASKRADLD l/D radius at central lambda
                // initialize the optical system
                PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 0);

                //     printf("-------------------------- STEP 0005a  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //            sleep(3);
                
                // computes or loads the piaa optics from the piaacmc structure
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                //   printf("-------------------------- STEP 0005b  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //          sleep(3);
                
                // initialize the optical system to be on axis
                PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);
                // printf("-------------------------- STEP 0005c  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
                //        sleep(3);

                // make the shapes again why?  *****************************************
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);

                /*               printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                                printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                                printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);

                                fflush(stdout);
                sleep(5);*/
                
                
                // focmMode controls which part of the FPM is propagated
                // if focmMode is a legal zone index, that zone index is propagated
                // here set focmMode beyond a legal zone index, so all zones are transparent and all
                // light including that which misses the FPM is propagated.
                // Later, we will subtract off the actual zone contributions, which will leave only
                // the light that misses the FPM.
                piaacmcsimul_var.focmMode = data.image[piaacmc[0].zonezID].md[0].size[0]+10;  // response for no focal plane mask
                optsyst[0].FOCMASKarray[0].mode = 1; // 1-fpm
                // compute the on-axis PSF to see what on-axis light goes around the FPM, return contrast in evaluation zone
                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                printf("val = %g\n", val);
                ID = image_ID("imvect");

                // FPMresp geometry:
                // first dimension (size[0]) is twice the number of evaluation points in the focal plane, giving Re and Im
                //      of the field at that evaluation point
                // second dimension (size[1]) is nbzones+1 zone indices, where nbzones is the number of mask zones (hexagons)
                // +1 because the first zone index stores the response to light that misses the FPM
                // WARNING: FPMresp size[1] is nbzones+1, as first vector stored is the response for light outside the mask
                // third dimension (size[2]) is wavelength


                // axis 0: eval pts (ii) - size = data.image[ID].md[0].size[0]
                // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
                // axis 3: lambda (k) - size = piaacmc[0].nblambda

                // allocate the combined FPMresp 3D array
                // ID is the "imvect" array created by PIAACMCsimul_computePSF and contains the pixels in the
                // evaluation set as a 1D vector (0th index) per wavelength
                IDfpmresp = create_3Dimage_ID_double("FPMresp", data.image[ID].md[0].size[0], piaacmc[0].focmNBzone+1, piaacmc[0].nblambda);
                //     list_image_ID();
                //    sleep(100);

                // light outside mask
                for(k=0; k<piaacmc[0].nblambda; k++) // loop over wavelengths
                    for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // loop over evaluation points
                        // set the 0th zone to be the light from the above on-axis PSF computation with a
                        // black FPM on the evaluation pixels in "imvect"
                        data.image[IDfpmresp].array.D[k*(piaacmc[0].focmNBzone+1)*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];


                if(piaacmcsimul_var.PIAACMC_FPMresp_thread > piaacmcsimul_var.PIAACMC_FPMresp_mp-1) // if we're the parent, combine files
                {
                    if((IDv=variable_ID("PID"))!=-1)
                        index = (long) data.variable[IDv].value.f+0.01;
                    // this file is looked for in the bash script, which waits for this
                    // file to spawn the tmux child processes
                    sprintf(command, "touch initOK_%ld", index);
                    printf("EXECUTING : %s\n", command);
                    r = system(command);
                    
                    // now the tmux children have been kicked off by the bash script and are
                    // computing their partial FPMresp files

                    // begin combining the partial files into the final FPMresp file when ready
                    printf("COMBINING FILES\n");
                    fflush(stdout);
                    sprintf(fnamet, "%s/FPMthreadstatus.txt", piaacmcsimul_var.piaacmcconfdir);
                    fpt = fopen(fnamet, "w");
                    fclose(fpt);

                    // we now wait for the children to complete their partial FPM resp files, then combine them
                    for(thr=0; thr<piaacmcsimul_var.PIAACMC_FPMresp_mp; thr++)
                    {
                        printf("thr = %ld\n", thr);
                        fflush(stdout);
                        
                        // each child creates an FPMresp...thread*.fits.tmp file, which is moved to
                        // FPMresp...thread*.fits (set as fname in the next sprintf) when the child is done
                        // signaling to the parent process that this part is ready to ingest.

                        ID1 = -1;
                        while(ID1 == -1) // wait for the partial FPMresp file from each child
                        {
                            // name of final child partial FPMresp file
                            sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d_mp%02ld_thread%02ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda, piaacmcsimul_var.PIAACMC_FPMresp_mp, thr);
                             
                            printf("Waiting for file \"%s\" ...\n", fname);
                            fflush(stdout);
                            // update thread status file
                            fpt = fopen(fnamet,"a");
                            fprintf(fpt, "Process %ld (thread %ld) --- Waiting for file \"%s\" ...\n", (long) getpid(), thr, fname);
                            fclose(fpt);
                            sleep(1.0);

                            // safely remove image with this name
                            delete_image_ID("tmpFPMresp");
                          //  list_image_ID();
                            // try to load the final child's partial FPMresp file
                            ID1 = load_fits(fname, "tmpFPMresp", 1);
                          //  list_image_ID();
                        }
                        // we found this child's partial FPMresp file!
                        // now insert it into our combined FPMresp file

                        fpt = fopen(fnamet,"a");
                        fprintf(fpt, "READING %s\n", fname);
                        fclose(fpt);

                        /*     list_image_ID();

                             printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                             printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                             printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);
                             printf("data.image[ID].md[0].size[0] =                   %ld\n", data.image[ID].md[0].size[0]);
                             printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n", data.image[piaacmc[0].zonezID].md[0].size[0]);
                             fflush(stdout);
                             sleep(100);
                        */
                        
                        assert( ID1 != -1 );
                        if(ID1!=-1) // this should always be true
                        {
                            mzstep = piaacmcsimul_var.PIAACMC_FPMresp_mp; // total number of tmux threads
                            mzoffset = thr; // the thread number of the child that just delivered its result
                            // insert the partial result of child thr into the combined FPMresp array
                            // be sure to skip the first line 'cause we already set it to be the light the went around the FPM
                            for(mz=1+mzoffset; mz<piaacmc[0].focmNBzone+1; mz+=mzstep) // loop over zone, do every PIAACMC_FPMresp_mp line
                            {

                                printf("mz = %ld    %ld %ld\n", mz, IDfpmresp, ID1);
                                fflush(stdout);
                                for(k=0; k<piaacmc[0].nblambda; k++) // for each wavelenth
                                    for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // for each evaluation point
                                    {
                                        // index of this evaluation point and wavelength and zone
                                        // tmpl1 = k*(nzones+1)*nEvaluationPoints) + zoneIndex*nEvaluationPoints + evaluationPoint
                                        tmpl1 = k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii;
                                        // set the combined array value from the partial file (both are same shape and size, of course)
                                        data.image[IDfpmresp].array.D[tmpl1] = data.image[ID1].array.D[tmpl1];
                                        // subtract the current zone value from the first zone line, which contained all light
                                        // (with no mask).  Eventually this will contain only light that misses the FPM.
                                        data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID1].array.D[tmpl1];
                                    }
                            }
                            delete_image_ID("tmpFPMresp"); // we're dont with the partial array, so delete it
                        }
                   
                    }
                    // write out the current state of the combined FPMresp file
                   sprintf(fname, "!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);


                    save_fits("FPMresp", fname);
                    // remove the child's .tmp file just in case (it should no longer exist 'cause we renamed, not copied, the .tmp file)
                    sprintf(command, "rm %s/FPMresp*.fits.tmp", piaacmcsimul_var.piaacmcconfdir);
                    r = system(command);
                }
                else // we're a tmux child, so compute the response for our portion
                {
                    // name of the child partial FPMresp file (to become .tmp)
                   sprintf(fname,"!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d_mp%02ld_thread%02ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda, piaacmcsimul_var.PIAACMC_FPMresp_mp, piaacmcsimul_var.PIAACMC_FPMresp_thread);


                    // diagnostic file to make sure the child is working with the right zones
                    sprintf(fnametmp, "!%s/fpmzmap_thread%02ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMresp_thread);
                    save_fits("fpmzmap", fnametmp);

                    printf("Making component %ld / %ld\n", piaacmcsimul_var.PIAACMC_FPMresp_thread, piaacmcsimul_var.PIAACMC_FPMresp_mp);
                    fflush(stdout);
                    piaacmcsimul_var.WRITE_OK = 0;
                    // for each FPM zone, compute the response
                    // skip the first one 'cause it is not computed by the children
                    for(mz=1+mzoffset; mz<piaacmc[0].focmNBzone+1; mz+=mzstep)
                    {
                        piaacmcsimul_var.focmMode = mz;  // focmMode can be a zone index, in which case operations are on that zone
                        optsyst[0].FOCMASKarray[0].mode = 0; // direct focal plane mask response

                        // default the reference to the 4th element
                        // but look for the element called "opaque mask at PIAA elem 1"
                        elem0 = 4;
                        for(elem=0; elem<optsyst[0].NBelem; elem++)
                        {
                            if(strcmp("opaque mask at PIAA elem 1", optsyst[0].name[elem])==0)
                            {
                                elem0 = elem;
                                printf("opaque mask at PIAA elem 1 = %ld\n", elem);
                            }
                             // raise an alarm if "opaque mask at PIAA elem 1" is not found *************************************
                        }

                        optsyst[0].keepMem[elem0] = 1; // keep it in memory

                        printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                        printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                        printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);
                        printf("data.image[ID].md[0].size[0] =                   %ld\n", (long) data.image[ID].md[0].size[0]);
                        printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n", (long) data.image[piaacmc[0].zonezID].md[0].size[0]);
                        fflush(stdout);

                        // compute the on-axis PSF
                        val = PIAACMCsimul_computePSF(0.0, 0.0, elem0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                        // The PSF result for the evaluation points is put in array "imvect" which previously was
                        // assigned to another PSF result.
                        // Should put another ID = image_ID("imvect") here *************************************

                        // set the response of this zone from the PSF result
                        for(k=0; k<piaacmc[0].nblambda; k++) // loop over wavelength
                            for(ii=0; ii<data.image[ID].md[0].size[0]; ii++) // loop over evaluation points
                            {
                                // see previous example for explanation of indexing
                                // save response, which is just the value of the on-axis PSF at each evaluation point
                                data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + mz*data.image[ID].md[0].size[0] + ii] = data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                                if(piaacmcsimul_var.PIAACMC_FPMresp_mp==1) // if we're single threaded (no children)
                                    // subtract the current zone value from the first zone line, which contained all light
                                    // (with no mask).  Eventually this will contain only light that misses the FPM.
                                    data.image[IDfpmresp].array.D[k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*data.image[ID].md[0].size[0] + ii] -= data.image[ID].array.F[k*data.image[ID].md[0].size[0]+ii];
                            }


                        printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname1);
                        fflush(stdout);
                        // fname1 is the .tmp name
                        save_fits("FPMresp", fname1);
                        printf("Done \n");
                        fflush(stdout);
                    }
                    
                    //*************************** make up our minds about single threading, in the meantime say we don't support it

                    // partial file complete!  move it to the final file name so parent can see it
                    printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname2);
                    fflush(stdout);
                    // fname2 is the final name
                    save_fits_atomic("FPMresp", fname2);
                    printf("Done \n");
                    fflush(stdout);
                    piaacmcsimul_var.WRITE_OK = 0;

                 //   sprintf(command, "mv %s %s", fname1, fname2);
                 //   r = system(command);
                }
            }
            else
                printf("File \"%s\" or \"%s\" exists\n", fname, fnamecomb);
        }
        piaacmcsimul_var.focmMode = -1;
	
	
	
	return 0;
}










    /** 
     * ---
     * 
     *  ## Mode 13: Optimize focal plane mask zones only
     * 
     * Uses "fast" mode:
     * 
     * After mode 11, we can use the (complex) light propagated from each zone to compute the impact of
     * any thickness (sag) of that zone: the zone thickness induces a phase rotation for that zone,
     * which is applied to the unobstructed light from that zone as a complex rotation.
     * 
     * The search is via steepest descent from random starting points.
     * 
     * This mode only sets up the optimization that actually happens after exiting the switch statement if piaacmcsimul_var.LINOPT = 1 (as does mode 40)
     * 
     */


int PIAACMCsimul_exec_mode13()
{
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;  	
	
    char stopfile[500];	
	
    int REGFPMSAG = 0; // regularization for FPM sag
 //   float fpmsagreg_coeff = 1.0e-8;
//    float fpmsagreg_coeff_alpha = 1.0;

    int paramtype[10000]; // FLOAT or DOUBLE
    double *paramval[10000]; // array of pointers, double
//    float *paramvalf[10000]; // array of pointers, float
//    double paramrefval[10000];

    double paramdelta[10000];
    double parammaxstep[10000]; // maximum single iteration step
    double parammin[10000]; // minimum value
    double parammax[10000]; // maximum value

long IDfpmresp;
long k;

	long IDstatus;
    char fname[1000]; 
    
	FILE *fp;
	
	long elem;
	int ret;
    double tmplf1, tmplf2;	
    int tmpd1;     
 
    long NBiter = 1000; 
	
	long ID;
	
	long mz;
	long number_param;
	
	        printf("=================================== mode 013 ===================================\n");


        // set the name of the stopfile
        sprintf(stopfile, "%s/stoploop13.txt", piaacmcsimul_var.piaacmcconfdir);



        // get cli variables

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        


        // FPM sag regularization control flag, do if == 1
        REGFPMSAG = 1; // default
        if((IDv=variable_ID("REGFPMSAG"))!=-1)
            REGFPMSAG = (long) data.variable[IDv].value.f+0.01;

     
  

        // set current state for statistical tracking
        data.image[IDstatus].array.UI16[0] = 0;

        // usual initialization
        PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);

        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);

        // set current state for statistical tracking
        data.image[IDstatus].array.UI16[0] = 1;

        // tracking diagnostic, giving the total flux in each plane
        sprintf(fname,"%s/flux.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "r");
        for(elem=0; elem<optsyst[0].NBelem; elem++)
        {
            ret = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
            optsyst[0].flux[elem] = tmplf1/tmpd1*optsyst[0].nblambda;   // scale flux to current number of lambda
        }
        fclose(fp);


        piaacmcsimul_var.LINOPT = 1; // perform linear optimization after the switch exits
        // get the number of iterations
        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            NBiter = (long) data.variable[IDv].value.f+0.01;
        else
            NBiter = 50; // default number of iterations

        // set current state for statistical tracking
        data.image[IDstatus].array.UI16[0] = 2;

        // get the FPMresp array computed in mode 11
        sprintf(fname, "%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);



        IDfpmresp = load_fits(fname, "FPMresp", 1);

        piaacmcsimul_var.vsize = data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2
        // make an array that holds the resulting light for evaluation point given the FPM solution, for each wavelenth
        ID = create_2Dimage_ID("imvect1", piaacmcsimul_var.vsize, piaacmc[0].nblambda);

        // allocate arrays for fast routine
        // define convenient array variables
        piaacmcsimul_var.fpmresp_array = data.image[IDfpmresp].array.D;
        piaacmcsimul_var.zonez_array = data.image[piaacmc[0].zonezID].array.D;
        // allocate derivative of phase against thickness array
        piaacmcsimul_var.dphadz_array = (double*) malloc(sizeof(double)*piaacmc[0].nblambda);
        // compute this derivative
        for(k=0; k<piaacmc[0].nblambda; k++)
        {
            // OPTICSMATERIALS_pha_lambda computes change in phase per unit thickness at specified wavelength
            // second arg is thickness, so 1.0 meters determines result per meter thickness
            piaacmcsimul_var.dphadz_array[k] = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, 1.0, optsyst[0].lambdaarray[k]);
            printf("%ld  %g %g\n", k, optsyst[0].lambdaarray[k], piaacmcsimul_var.dphadz_array[k]);
        }
        piaacmcsimul_var.outtmp_array = (double*) malloc(sizeof(double)*(piaacmcsimul_var.vsize*piaacmc[0].nblambda+data.image[piaacmc[0].zonezID].md[0].size[0]));

        // do the fast optimization using the results of mode 11
        piaacmcsimul_var.computePSF_FAST_FPMresp = 1;

        // set current state for statistical tracking
        data.image[IDstatus].array.UI16[0] = 3;

        // read the contrast normalization factor into CnormFactor
        sprintf(fname, "%s/CnormFactor.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "r");
        ret = fscanf(fp, "%lf", &piaacmcsimul_var.CnormFactor);
        fclose(fp);
        // for each zone, add a random offset in range +- MODampl
        // this randomizes the starting point for each zone
        // data.image[piaacmc[0].zonezID].array.D[k] is set in PIAACMCsimul_run()
        for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
            data.image[piaacmc[0].zonezID].array.D[k] += piaacmcsimul_var.MODampl*(1.0-2.0*ran1());

        // set up optimization parameters for each zone
        // uses abstract specification of optimization parameters called paramval, ...
        number_param = 0;
        for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
        {
            // parameter type
            paramtype[number_param] = _DATATYPE_DOUBLE;
            // value: sag of each zone
            paramval[number_param] = &data.image[piaacmc[0].zonezID].array.D[mz];
            // derivative step size
            paramdelta[number_param] = 3.0e-9;
            // max parameter step size
            parammaxstep[number_param] = 1.0e-6;
            // max and min allowed values for parameter
            parammin[number_param] = piaacmc[0].fpmminsag;
            parammax[number_param] = piaacmc[0].fpmmaxsag;
            // move on to next parameter
            number_param++;
        }
        piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES = 1; // for fast execution using analytic derivatives

        // set current state for statistical tracking
        data.image[IDstatus].array.UI16[0] = 4;
        // Now on to the actual optimization, after exit from the switch statement
        // I hope you have a lot of time...
	
	return 0;
}





















    /**
     * ---
     * 
     * ## Mode 40: Optimize PIAA optics shapes (and focal plane mask transmission for idealized PIAACMC)
     * 
     * 
     */ 
     
int PIAACMCsimul_exec_mode40()
{    
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;  	
    long NBiter = 1000;	
	
	long number_param = 0;
	
    int paramtype[10000]; // FLOAT or DOUBLE
    double *paramval[10000]; // array of pointers, double
    float *paramvalf[10000]; // array of pointers, float
//    double paramrefval[10000];

    double paramdelta[10000];
    double parammaxstep[10000]; // maximum single iteration step
    double parammin[10000]; // minimum value
    double parammax[10000]; // maximum value
    
    long mz, k;
    long ID_CPAfreq;
    long kmaxC, kmaxF;
    
        // OPTIMIZATION PARAMETERS
    int REGPIAASHAPES = 0;

    float piaa0C_regcoeff = 0.0e-7; // regularization coeff
    float piaa1C_regcoeff = 0.0e-7; // regularization coeff
    float piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
    float piaa1C_regcoeff_alpha = 1.0; // regularization coeff power

    float piaa0F_regcoeff = 0.0e-7; // regularization coeff
    float piaa1F_regcoeff = 0.0e-7; // regularization coeff
    float piaa0F_regcoeff_alpha = 1.0; // regularization coeff power
    float piaa1F_regcoeff_alpha = 1.0; // regularization coeff power

    int REGFPMSAG = 0; // regularization for FPM sag
 //   float fpmsagreg_coeff = 1.0e-8;
//    float fpmsagreg_coeff_alpha = 1.0;

    	
	
        printf("=================================== mode 040 ===================================\n");
        //		piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
 
     if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        

 
 
        piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);
        printf("PIAACMC_fpmtype = %d\n", piaacmcsimul_var.PIAACMC_fpmtype);


        piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
        PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);
        //       printf("data.image[piaacmc[0].zoneaID].array.D[0] = %lf\n", data.image[piaacmc[0].zoneaID].array.D[0]);
        //        sleep(10);



        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


        if(0) // TEST
        {
			double valref;
			
            valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
            printf("valref = %g\n", valref);
            printf("EXEC CASE 40 COMPLETED\n");
        }
        else
        {
			piaacmcsimul_var.LINOPT = 1; // perform linear optimization
			
            if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
                NBiter = (long) data.variable[IDv].value.f+0.01;
            else
                NBiter = 1000;


            kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
                kmaxC = (long) data.variable[IDv].value.f+0.01;
            if(kmaxC>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
                kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];


            kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
            if((IDv=variable_ID("PIAACMC_maxoptFterm"))!=-1)
                kmaxF = (long) data.variable[IDv].value.f+0.01;
            if(kmaxF>data.image[piaacmc[0].piaa0FmodesID].md[0].size[0])
                kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];




            // PIAA shapes regularization

            REGPIAASHAPES = 0; // default
            if((IDv=variable_ID("REGPIAASHAPES"))!=-1)
                REGPIAASHAPES = (long) data.variable[IDv].value.f+0.01;

            piaa0C_regcoeff = 0.0e-7; // regularization coeff
            piaa1C_regcoeff = 0.0e-7; // regularization coeff
            if((IDv=variable_ID("REGPIAA_C_COEFF"))!=-1)
            {
                piaa0C_regcoeff = data.variable[IDv].value.f;
                piaa1C_regcoeff = data.variable[IDv].value.f;
            }

            piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
            piaa1C_regcoeff_alpha = 1.0; // regularization coeff power
            if((IDv=variable_ID("REGPIAA_C_ALPHA"))!=-1)
            {
                piaa0C_regcoeff_alpha = data.variable[IDv].value.f;
                piaa1C_regcoeff_alpha = data.variable[IDv].value.f;
            }

            piaa0F_regcoeff = 0.0e-7; // regularization coeff
            piaa1F_regcoeff = 0.0e-7; // regularization coeff
            if((IDv=variable_ID("REGPIAA_F_COEFF"))!=-1)
            {
                piaa0F_regcoeff = data.variable[IDv].value.f;
                piaa1F_regcoeff = data.variable[IDv].value.f;
            }

            piaa0F_regcoeff_alpha = 1.0; // regularization coeff power
            piaa1F_regcoeff_alpha = 1.0; // regularization coeff power
            if((IDv=variable_ID("REGPIAA_F_ALPHA"))!=-1)
            {
                piaa0F_regcoeff_alpha = data.variable[IDv].value.f;
                piaa1F_regcoeff_alpha = data.variable[IDv].value.f;
            }

            if(REGPIAASHAPES==1)
            {
				printf("loading CPA modes frequ\n");
				fflush(stdout);
                ID_CPAfreq = image_ID("cpamodesfreq");
                if(ID_CPAfreq == -1)
                    ID_CPAfreq = load_fits("cpamodesfreq.fits", "cpamodesfreq", 2);
            }
			

            // FPM SAG regularization

            REGFPMSAG = 1; // default
            if((IDv=variable_ID("REGFPMSAG"))!=-1)
                REGFPMSAG = (long) data.variable[IDv].value.f+0.01;




            number_param = 0;

            if(piaacmcsimul_var.PIAACMC_fpmtype==0) // ideal mask
            {
                paramtype[number_param] = _DATATYPE_DOUBLE;
                paramval[number_param] = &data.image[piaacmc[0].zoneaID].array.D[0];
                paramdelta[number_param] = 1.0e-3;
                parammaxstep[number_param] = 2.0e-1;
                parammin[number_param] = -1.0;
                parammax[number_param] = 1.0;
                number_param++;
            }
            else // real physical mask
            {
                if(variable_ID("PIAACMC_mzOPT")!=-1) // optimize zones
                {
                    for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
                    {
                        paramtype[number_param] = _DATATYPE_DOUBLE;
                        paramval[number_param] = &data.image[piaacmc[0].zonezID].array.D[mz];
                        paramdelta[number_param] = 1.0e-9;
                        parammaxstep[number_param] = 5.0e-8;
                        parammin[number_param] = -2.0e-6;
                        parammax[number_param] = 2.0e-6;
                        number_param++;
                    }
                }
            }


            for(k=0; k<kmaxC; k++)
            {
                paramtype[number_param] = _DATATYPE_FLOAT;
                paramvalf[number_param] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
                paramdelta[number_param] = 1.0e-10;
                parammaxstep[number_param] = 1.0e-7;
                parammin[number_param] = -1.0e-3;
                parammax[number_param] = 1.0e-3;
                number_param++;
            }

            for(k=0; k<kmaxC; k++)
            {
                paramtype[number_param] = _DATATYPE_FLOAT;
                paramvalf[number_param] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
                paramdelta[number_param] = 1.0e-10;
                parammaxstep[number_param] = 1.0e-7;
                parammin[number_param] = -1.0e-3;
                parammax[number_param] = 1.0e-3;
                number_param++;
            }

            for(k=0; k<kmaxF; k++)
            {
                paramtype[number_param] = _DATATYPE_FLOAT;
                paramvalf[number_param] = &data.image[piaacmc[0].piaa0FmodesID].array.F[k];
                paramdelta[number_param] = 1.0e-10;
                parammaxstep[number_param] = 1.0e-7;
                parammin[number_param] = -1.0e-3;
                parammax[number_param] = 1.0e-3;
                number_param++;
            }

            for(k=0; k<kmaxF; k++)
            {
                paramtype[number_param] = _DATATYPE_FLOAT;
                paramvalf[number_param] = &data.image[piaacmc[0].piaa1FmodesID].array.F[k];
                paramdelta[number_param] = 1.0e-10;
                parammaxstep[number_param] = 1.0e-7;
                parammin[number_param] = -1.0e-3;
                parammax[number_param] = 1.0e-3;
                number_param++;
            }

            piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 1;
            piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 1;
        }

	return 0;
}










int PIAACMCsimul_exec_mode100()
{
	long IDopderrC;
	long nbOPDerr;
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;  	
    double ldoffset;	
	double valref;
	
	long ID;
	char fname[1000];
	long xsize, ysize, zsize;
	
    float *peakarray;
    float avpeak;

	long kk, ii, jj;
	double val;

    double focscale;
    double xc, yc, rc;
	
    double *eval_contrastCurve;
    double *eval_COHcontrastCurve;
    double *eval_INCcontrastCurve;
    double *eval_contrastCurve_cnt;
    double eval_sepstepld = 0.2; // in l/D
    double eval_sepmaxld = 20.0; // in l/D
    long eval_sepNBpt;
	
	long ri;
	FILE *fp;
	long IDps, IDps_re, IDps_im, IDps_COH, IDps_INC, IDre, IDim, IDopderr, IDrc;
	
		float *rcarray;
	float drc;
	long rcbin;
	long rcbinmax;
	float rcmin, rcmax;
	long rc_cnt;
	float p50, p20;
	
	long NBpt;
	
    double aveC;
    double aveC_COH;
    double aveC_INC;
    long aveCcnt;	
	
	long OPDmode;
	long size;
	
        printf("=================================== mode 100 ===================================\n");
		
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        

		

		// measure sensitivity to errors
	//	printf("Loading (optional) OPDerr file\n");
	//	fflush(stdout);
		// load an error if it exists
		IDopderrC = image_ID("OPDerrC");
		if(IDopderrC == -1)
			IDopderrC = load_fits("OPDerrC.fits", "OPDerrC", 0);
	

		if(IDopderrC != -1)
		{
			nbOPDerr = data.image[IDopderrC].md[0].size[2];  // number of error arrays
			//printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
			//fflush(stdout);
		}
		else
		{
			//printf("NO OPD ERROR MODES\n");
			//fflush(stdout);
			nbOPDerr = 0;
		}

	
  
		printf("Will add optional OPD error modes (%ld modes)\n", nbOPDerr);
		fflush(stdout);
		


        piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);


	

        piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
        PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);





        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        //        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0);
        //       printf("valref = %g\n", valref);

        //exit(0);




        ldoffset = 0.01; // default
        if((IDv=variable_ID("PIAACMC_ldoffset"))!=-1)
            ldoffset = data.variable[IDv].value.f;

        printf("ldoffset = %f\n", ldoffset);

		// compute off-axis POINT source
        valref = PIAACMCsimul_computePSF(5.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_x50_y00.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        //load_fits(fname, "psfi");




        ID = image_ID("psfi0");
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];
        zsize = data.image[ID].md[0].size[2];
        peakarray = (float*) malloc(sizeof(float)*zsize);
        for(kk=0; kk<zsize; kk++)
        {
            peakarray[kk] = 0.0;
            for(ii=0; ii<xsize*ysize; ii++)
            {
                val = data.image[ID].array.F[kk*xsize*ysize+ii];
                if(val>peakarray[kk])
                    peakarray[kk] = val;
            }
        }
        avpeak = 0.0;
        for(kk=0; kk<zsize; kk++)
        {
            printf("peak %02ld  %10lf\n", kk, peakarray[kk]);
            avpeak += peakarray[kk];
        }
        avpeak /= zsize;
        free(peakarray);
        delete_image_ID("psfi0");


		// compute on-axis POINT source
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
        sprintf(fname,"!%s/psfi0_x00_y00.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        //load_fits(fname, "psfi");


        ID = image_ID("psfi0");




        /// compute contrast curve
        focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
        printf("focscale = %f\n", focscale);
        eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
        eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        aveC = 0.0;
        aveCcnt = 0;

        for(kk=0; kk<zsize; kk++)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    xc = 1.0*ii-0.5*xsize;
                    yc = 1.0*jj-0.5*ysize;
                    xc *= focscale;
                    yc *= focscale;
                    rc = sqrt(xc*xc+yc*yc);
                    ri = (long) (rc/eval_sepstepld-0.5);
                    if(ri<0)
                        ri = 0;
                    if(ri<eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if((rc>2.0)&&(rc<6.0))
                    {
                        aveC += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveCcnt++;
                    }
                }



        sprintf(fname, "%s/ContrastCurve_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d_tt000.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        fp = fopen(fname, "w");
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
        }
        fclose(fp);
        free(eval_contrastCurve);
        free(eval_contrastCurve_cnt);

       sprintf(fname, "%s/ContrastVal_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d_tt000.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


        fp = fopen(fname, "w");
        fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode);
        fclose(fp);




        // measure pointing sensitivity
        IDps = create_3Dimage_ID("starim", piaacmc[0].size, piaacmc[0].size, zsize);

		IDps_re = create_3Dimage_ID("starim_re", piaacmc[0].size, piaacmc[0].size, zsize);
		IDps_im = create_3Dimage_ID("starim_im", piaacmc[0].size, piaacmc[0].size, zsize);		
		
		IDps_COH = create_3Dimage_ID("starimCOH", piaacmc[0].size, piaacmc[0].size, zsize);
		IDps_INC = create_3Dimage_ID("starimINC", piaacmc[0].size, piaacmc[0].size, zsize);
		
		NBpt = 0;
		
        valref = 0.25*PIAACMCsimul_computePSF(ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_p0.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        for(ii=0; ii<xsize*ysize*zsize; ii++){
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
			data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
			data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
		}
		NBpt++;
		delete_image_ID("psfre0");
		delete_image_ID("psfim0");
		
        valref += 0.25*PIAACMCsimul_computePSF(-ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_m0.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        for(ii=0; ii<xsize*ysize*zsize; ii++){
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
			data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
			data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
		}
		NBpt++;
		delete_image_ID("psfre0");
		delete_image_ID("psfim0");

        valref += 0.25*PIAACMCsimul_computePSF(0.0, ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_0p.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        for(ii=0; ii<xsize*ysize*zsize; ii++){
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
			data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
			data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
		}
		NBpt++;
		delete_image_ID("psfre0");
		delete_image_ID("psfim0");

        valref += 0.25*PIAACMCsimul_computePSF(0.0, -ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname,"!%s/psfi0_0m.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);
        ID = image_ID("psfi0");
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        for(ii=0; ii<xsize*ysize*zsize; ii++){
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
			data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
			data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
		}
		NBpt++;
		delete_image_ID("psfre0");
		delete_image_ID("psfim0");

      
	
	
	 // measure sensitivity to errors

//	printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
//	fflush(stdout);
	// add error modes if any
		for(OPDmode=0; OPDmode < nbOPDerr; OPDmode++)
		{
			size = data.image[IDopderrC].md[0].size[0];
			IDopderr = create_2Dimage_ID("opderr", size, size);
            // "opderr" is a standard name read by PIAACMCsimul_init
			for(ii=0;ii<size*size;ii++)
				data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size*size*OPDmode + ii];			
			PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // add error to the data
			PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
			sprintf(fname, "!%s/psfi0_opderr%02ld.fits", piaacmcsimul_var.piaacmcconfdir, OPDmode);
			save_fits("psfi0", fname);
			delete_image_ID("opderr");
			ID = image_ID("psfi0");
			IDre = image_ID("psfre0");
			IDim = image_ID("psfim0");
			for(ii=0; ii<xsize*ysize*zsize; ii++){
				data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
				data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
				data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
			}
			NBpt++;
			delete_image_ID("psfre0");
			delete_image_ID("psfim0");
		}
		
		
		for(ii=0; ii<xsize*ysize*zsize; ii++){
            data.image[IDps].array.F[ii] /= NBpt;
			data.image[IDps_re].array.F[ii] /= NBpt; // average re
			data.image[IDps_im].array.F[ii] /= NBpt; // average im
		}
		
		
		for(ii=0; ii<xsize*ysize*zsize; ii++){
			data.image[IDps_COH].array.F[ii] = data.image[IDps_re].array.F[ii]*data.image[IDps_re].array.F[ii] + data.image[IDps_im].array.F[ii]*data.image[IDps_im].array.F[ii];
			data.image[IDps_INC].array.F[ii] = data.image[IDps].array.F[ii] - data.image[IDps_COH].array.F[ii];
		}



        sprintf(fname, "!%s/psfi0_starim.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("starim", fname);




        sprintf(fname, "!%s/psfi0_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
        save_fits("starim", fname);

        sprintf(fname, "!%s/psfi0INC_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
        save_fits("starimINC", fname);

        sprintf(fname, "!%s/psfi0COH_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
        save_fits("starimCOH", fname);




        /// compute contrast curve
        /// measure average contrast value, 2-6 lambda/D
        focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
        printf("focscale = %f\n", focscale);
        eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
        eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_COHcontrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_INCcontrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
        eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_COHcontrastCurve[ri] = 0.0;
            eval_INCcontrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        aveC = 0.0;
        aveC_INC = 0.0;
        aveC_COH = 0.0;
        aveCcnt = 0;


		IDrc = create_3Dimage_ID("_tmp_rc", xsize, ysize, zsize);
		rcarray = (float*) malloc(sizeof(float)*xsize*ysize*zsize);
		

        for(kk=0; kk<zsize; kk++)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    xc = 1.0*ii-0.5*xsize;
                    yc = 1.0*jj-0.5*ysize;
                    xc *= focscale;
                    yc *= focscale;
                    rc = sqrt(xc*xc+yc*yc);
                    data.image[IDrc].array.F[kk*xsize*ysize+jj*xsize+ii] = rc;
                    
                    ri = (long) (rc/eval_sepstepld-0.5);
                    if(ri<0)
                        ri = 0;
                    if(ri<eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_COHcontrastCurve[ri] += data.image[IDps_COH].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_INCcontrastCurve[ri] += data.image[IDps_INC].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if((rc>2.0)&&(rc<6.0))
                    {
                        aveC += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveC_COH += data.image[IDps_COH].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveC_INC += data.image[IDps_INC].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        aveCcnt++;
                    }
                }

		sprintf(fname, "%s/ContrastCurveP_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

		fp = fopen(fname, "w");
		drc = 0.2;
		rcbinmax = 100;
		for(rcbin=0;rcbin<rcbinmax;rcbin++)
			{
				rcmin = drc * rcbin;
				rcmax = drc * (rcbin+1);
				
				rc_cnt = 0;
				
				for(kk=0; kk<zsize; kk++)
				for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
				{
				rc = data.image[IDrc].array.F[kk*xsize*ysize+jj*xsize+ii];
				if((rc>rcmin)&&(rc<rcmax))
					{
						rcarray[rc_cnt] = data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
						rc_cnt++;
					}
				}
				quick_sort_float(rcarray, rc_cnt);
			
				p50 = rcarray[(long) (0.5*rc_cnt)];
				p20 = rcarray[(long) (0.2*rc_cnt)];
				fprintf(fp, "%5.2f %12g %12g\n", drc*(rcbin+0.5), p50, p20);
			}


		free(rcarray);
		delete_image_ID("_tmp_rc");
		fclose(fp);
		
		delete_image_ID("starimINC");
		delete_image_ID("starimCOH");


        sprintf(fname, "%s/ContrastCurve_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);


        fp = fopen(fname, "w");
        for(ri=0; ri<eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            eval_COHcontrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            eval_INCcontrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
            fprintf(fp, "%10f %10g %10g %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri], eval_COHcontrastCurve[ri], eval_INCcontrastCurve[ri]);
        }
        fclose(fp);

        free(eval_contrastCurve);
        free(eval_COHcontrastCurve);
        free(eval_INCcontrastCurve);
        free(eval_contrastCurve_cnt);

        sprintf(fname, "%s/ContrastVal_extsrc%2ld_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

        fp = fopen(fname, "w");
        fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %7.3f %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) (1000.0*ldoffset), piaacmc[0].fpmsagreg_coeff, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode);
        fclose(fp);

return 0;
}















    /**
     * ---
     * 
     * ## Mode 101: Measure transmission as a function of angular separation
     * 
     */ 
int PIAACMCsimul_exec_mode101()
{
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;  	
     double ldoffset;
	double valref;
	char fname[1000];
	char fnametransm[1000];

        printf("=================================== mode 101 ===================================\n");

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }        


        printf("101: transm as a function of angular separation  ldoffset = %f\n", ldoffset);

        piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
            piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);

        piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
        PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);

        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
        sprintf(fname,"!%s/psfi0test_x00_y00.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("psfi0", fname);


       sprintf(fnametransm, "%s/transmCurve_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

		FILE *fpt;
        fpt = fopen(fnametransm, "w");
        fclose(fpt);
		
		double stepld = 0.001;
        double xld;
        for(xld=0.0; xld<10.0; xld+=stepld)
        {
			long ID;
			long xsize, ysize, zsize;
			double val;
			
            valref = PIAACMCsimul_computePSF(xld, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
            ID = image_ID("psfi0");
            //            sprintf(fname, "!psfi0transm_%04.1f.fits", xld);
            //           save_fits("psfi0", fname);
            printf("ID = %ld\n", ID);
            xsize = data.image[ID].md[0].size[0];
            ysize = data.image[ID].md[0].size[1];
            zsize = data.image[ID].md[0].size[2];
            printf("image size = %ld %ld %ld\n", xsize, ysize, zsize);
            val = 0.0;
            
            long kk;
            for(kk=0; kk<zsize; kk++)
            {
				long ii, jj;
                for(ii=0; ii<xsize; ii++)
                    for(jj=0; jj<ysize; jj++)
                    {
						double dx, dy;
						
                        dx = 1.0*ii-0.5*xsize;
                        dy = 1.0*jj-0.5*ysize;
                        if((dx*dx+dy*dy)<30.0*30.0)
                            val += data.image[ID].array.F[kk*xsize*ysize+jj*ysize+ii];
                    }
            }
            val /= zsize;

            fpt = fopen(fnametransm, "a");
            fprintf(fpt, "%10f %.18f\n", xld, val);
            fclose(fpt);
            delete_image_ID("psfi0");

            stepld = 0.001;
            stepld += 0.1*xld;
            if(stepld>0.2)
                stepld = 0.2;
        }

	return 0;
}





















/**
 *
 * @brief Main simulation routine
 *
 * 
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 *
 */

int PIAACMCsimul_exec(const char *confindex, long mode)
{
    FILE *fp;
    char command[1000];
    char dirname[500];
    int paramtype[10000];
    double *paramval[10000]; 
    float *paramvalf[10000]; 

    double paramdelta[10000];
    double parammaxstep[10000]; 
    double parammin[10000]; 
    double parammax[10000]; 

    double paramdeltaval[10000];


    double valref, val0;
   long i, jj, ii;

    long IDv, ID, IDref;
    long IDmodes, IDmodes2D;
    long xsize = 0;
    long ysize = 0;
    long k;

    long iter;
    long NBiter = 1000;


    char fname[800];
    long IDm, ID1D, ID1Dref;
    long size1Dvec, size1Dvec0;




  
    int REGPIAASHAPES = 0;

    float piaa0C_regcoeff = 0.0e-7; 
    float piaa1C_regcoeff = 0.0e-7; 
    float piaa0C_regcoeff_alpha = 1.0; 
    float piaa1C_regcoeff_alpha = 1.0; 

    float piaa0F_regcoeff = 0.0e-7; 
    float piaa1F_regcoeff = 0.0e-7; 
    float piaa0F_regcoeff_alpha = 1.0; 
    float piaa1F_regcoeff_alpha = 1.0;


    int REGFPMSAG = 0; 

    int r;
    double val;

    double fpmradld = 0.95;  
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    long mz;

    long ii1, ii2, ki, kv, ki1;

    double scangain;
    double scanstepgain = 0.001;
    int linscanOK;
    double valold, oldval;
    double bestgain;
    double linoptgainarray[100];
    double linoptvalarray[100];
    int linoptlimflagarray[100];
    long NBlinoptgain;
    long kmax;


    double valtest;
    int ret;
    int kmaxC, kmaxF;
	long number_param = 0;

    long st;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function fpmeval_func;
    int status;

    // simulated annealing loop(s)
    long NBITER_SA = 1000000;
    long ITERMAX;
    long iter1;
    double SAcoeff;
    double amp1;
    double val1;
    int SA_MV = 0;
    double tmp1;
    double bestvalue;
    FILE *fpbest;

    // downhill simplex
    long NBITER_DS = 1000;
    long iterMax;
    int OK;
    int KEEP;
    double KEEPlimit = 1.0;
    double KEEPlimit1 = 1.0;
    double eps1;
    double mmsize;

    long elem;


    int iterOK;
    char stopfile[500];

    double alphareg;
    double bestval = 0.0;
    long IDoptvec = -1;

    double acoeff0, acoeff1, acoeff2;
    float scangainfact;
    int alphascaninit = 0;

    uint32_t *sizearray;
    long IDstatus;

    double valContrast;
    double tmp;
    int initbestval = 0;

	long ID_CPAfreq;





	#ifdef PIAASIMUL_LOGFUNC0
		sprintf(flogcomment, "%s %ld", confindex, mode);
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, flogcomment);
	#endif


    // Create status shared variable
    // this allows realtime monitoring of the code by other processes
    // sets status at different points in the code
    IDstatus = image_ID("stat_PIAACMCsimulexec");
    if(IDstatus == -1)
        IDstatus = read_sharedmem_image("stat_PIAACMCsimulexec");

    sizearray = (uint32_t*) malloc(sizeof(uint32_t)*2);
    sizearray[0] = 1;
    sizearray[1] = 1;
    IDstatus = create_image_ID("stat_PIAACMCsimulexec", 2, sizearray, _DATATYPE_UINT16, 1, 0);
    free(sizearray);



    //piaacmc = NULL; // set the pointer to the piaacmc structure to null

    // if the optical system pointer is empty, create an empty version
    if(optsyst==NULL)
    {
        optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));
        optsyst[0].SAVE = 1;
    }

    for(elem=0; elem<100; elem++)
        optsyst[0].keepMem[elem] = 0; // flag that says save this element for reuse


    // set the result directories
    sprintf( piaacmcsimul_var.piaacmcconfdir, "%s", confindex);
    sprintf( data.SAVEDIR, "%s", piaacmcsimul_var.piaacmcconfdir);






    optsyst[0].SAVE = piaacmcsimul_var.PIAACMC_save;

    // get variables from command line, possibly sets globals
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    // start a log of code mode entry/exit times
    sprintf(command, "echo \"%03ld     $(date)\" >> ./log/PIAACMC_mode_log.txt", mode);
    ret = system(command);
    printf("command = %s\n", command);



    // set the name of the stopfile
    sprintf(stopfile, "%s/stopmode%ld.txt", piaacmcsimul_var.piaacmcconfdir, mode);

    switch (mode) {
		
    case 0 :  
		PIAACMCsimul_exec_compute_image();
        printf("EXEC CASE 0 COMPLETED\n");
        fflush(stdout);
        break;

    case 1 : 
		PIAACMCsimul_exec_optimize_lyot_stop_position();
        break;

    case 2 : 
		PIAACMCsimul_exec_optimize_fpmtransmission();
        break;

    case 3 : 
		val = PIAACMCsimul_exec_mode03();
        break;

    case 4 :
        PIAACMCsimul_exec_mode04();
        break;

    case 5 :
		PIAACMCsimul_exec_mode05();
        break;

    case 11 :
		PIAACMCsimul_exec_mode11();
        break;

    case 13 : 
		PIAACMCsimul_exec_mode13();
        break;

    case 40 : 
		PIAACMCsimul_exec_mode40();
        break;

    case 100 : // evaluate current design: polychromatic contrast, pointing sensitivity
		PIAACMCsimul_exec_mode100();
        break;

    case 101 : 
		PIAACMCsimul_exec_mode101();
        break;


    case 300 : 
    /** 
     * ---
     * 
     * ## Mode 300: Import FPM configuration setting from parent directory
     * 
     */ 
        printf("=================================== mode 300 ===================================\n");

        /*  sprintf(command, "cp conf_MASKRADLD.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_FPMsectors.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_NBrings.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_nblambda.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_resolved.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);
          sprintf(command, "cp conf_extmode.txt %s/", piaacmcsimul_var.piaacmcconfdir);
          r = system(command);*/
        break;

    case 301 : // remove configuration settings
        printf("=================================== mode 301 ===================================\n");

        /*   sprintf(command, "mv %s/conf_MASKRADLD.txt %s/saveconf/conf_MASKRADLD.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_FPMsectors.txt %s/saveconf/conf_FPMsectors.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_NBrings.txt %s/saveconf/conf_NBrings.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_nblambda.txt %s/saveconf/conf_nblambda.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_resolved.txt %s/saveconf/conf_resolved.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);
           sprintf(command, "mv %s/conf_extmode.txt %s/saveconf/conf_extmode.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
           r = system(command);*/
        break;

    case 302 : // restore configuration settings
        printf("=================================== mode 302 ===================================\n");

        sprintf(command, "cp %s/saveconf/conf_*.txt %s/", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        r = system(command);
        break;


    default :
        printERROR(__FILE__,__func__,__LINE__, "mode not recognized");
        break;
    }



































    // linear optimization set up in modes 13 and 40
    //
    // the output parameters start as the evaluation zone values in "imvect"
    // If we are regularizing, we supplement the output parameters by adding
    // penalty terms as additional output parameters
    if(piaacmcsimul_var.LINOPT == 1) // linear optimization
    {
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 5;

        // Compute Reference on-axis performance contrast (valref)
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);


        // save the current configuration to the _linopt directory
        sprintf(dirname, "%s_linopt", piaacmcsimul_var.piaacmcconfdir);
        PIAACMCsimul_savepiaacmcconf(dirname);

        // import configuration from _linopt directory
        sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcsimul_var.piaacmcconfdir);
        r = system(command);

        // "cp -n" will only copy the file if destination does not exist
        // this will ensure that the .ref.fits files are the PIAA shapes before any linear optimization
        // these will be the "reference" shapes used to regularize PIAA shapes: the deviation from this reference will be kept small by the regularization

        sprintf(command, "cp -n %s/piaa0Cmodes.fits %s/piaa0Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa0Fmodes.fits %s/piaa0Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaa1Cmodes.fits %s/piaa1Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa1Fmodes.fits %s/piaa1Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaacmcparams.conf %s/piaacmcparams.ref.conf", dirname, dirname);
        r = system(command);


        // load PIAA reference shapes (used for regularization)
        sprintf(fname, "%s/piaa0Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Cmref", 1);

        sprintf(fname, "%s/piaa1Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Cmref", 1);

        sprintf(fname, "%s/piaa0Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Fmref", 1);

        sprintf(fname, "%s/piaa1Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Fmref", 1);

        // we have now saved the starting point of the optimization for future comparison
        // in the <piaacmcconfdir>_linopt directory

        // here we compute regularization value of piaashapes and store it in the val0 variable
        // if regularization is of PIAA shapes is ON, then val0 will be computed and added to the overal performance metric valref
        val0 = 1.0;

        // regularize the piaashapes via a penalty added to the reference contrast valref
        // The optimization minimizes the summed contrast + val0 + val1.
        // Regularization is via adding a constant val0 + val1 to the contrast we're minimizing
        // note that here we're setting output parameters.
        if(REGPIAASHAPES==1)
        {
            // first we compute the starting regularization constant
            val0 = 0.0;
            // index of the PIAA element 0 (first mirror) shapes via cosine modes
            ID = piaacmc[0].piaa0CmodesID;
            // index of PIAA shapes reference image
            IDref = image_ID("piaa0Cmref");
            if(IDref==-1)
            {   // error message if we get here?  ***************************************
                // if the reference image doesn't exist, create it
                IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                // initialize to zero shape
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            // This section of code does not actually fill in the regularization terms in the output vector
            // filling in is done later.  Here we are only computing the initial reference scalar objective value
            
            // for each cosine mode set the optimization parameter = cosine mode modified by regularization
            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            {
                // compute square of C*(deviation from reference)*(mode index)^(alpha)
                // so higher-index zones (higher spatial frequency) are more
                // heavily penalized
                tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                val0 += tmp*tmp;
            }


            // do the same for PIAA element 1
            ID = piaacmc[0].piaa1CmodesID;
            IDref = image_ID("piaa1Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            // get spatial frequency of each mode in cycles/aperture
            ID_CPAfreq = image_ID("cpamodesfreq");

            // do the same for PIAA element 0 and 1 for the Fourier modes
            // this time use the actual spatial frequency rather than mode index as proxy for frequency
            ID = piaacmc[0].piaa0FmodesID;
            IDref = image_ID("piaa0Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            ID = piaacmc[0].piaa1FmodesID;
            IDref = image_ID("piaa1Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
            {
                tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                val0 += tmp*tmp;
            }

            printf("VALREF = %g + %g -> %g\n", valref, val0, valref+val0);
            valref += val0;
        }


        // val1 is the regularization value for the focal plane mask sag values
        // same as above, but not dependent on position
        val1 = 1.0;
        if(REGFPMSAG == 1)
        {
            ID = piaacmc[0].zonezID;
            val1 = 0.0;
            for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
            {
                // compute the square of (sag/coeff)^alpha
                tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                val1 += tmp*tmp;
            }
            valref += val1;
        }
        // At this point all we've done is compute the overall performance metric including
        // regularization in valref.


        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 6;
        printf("================================ Reference = %g\n", valref);

        // copy imvect to vecDHref "vector dark hole reference"
        // vecDHref is the dark hole complex amplitude state at the beginning of the linear optimization
        // (this dark hole is nominally a full annulus from 1.5 to ~8 lambda/D, created at the top
        // of PIAACMCsimul_computePSF with size controlled by scoringIWA and scoringOWA
        // the corresponding performance metric is valref
        chname_image_ID("imvect", "vecDHref"); // note: imvect was computed by PIAACMCsimul_computePSF called ~150 lines above
        ID = image_ID("vecDHref"); // ID changed identity
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];



        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 7;

        // now we will just determine the size of the size of the
        // optimization vectors that we will actually fill in later
        // save vecDHref initial state as a reference
        sprintf(fname, "!%s/vecDMref.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("vecDHref", fname);
        // get and copy size of vecDHref, 'cause we're manipulating size1Dvec
        size1Dvec = data.image[ID].md[0].nelement;
        size1Dvec0 = size1Dvec;

        // PIAA shapes regularization
        // if regularization is turned on, the size of the evaluation vector is increased to include PIAA shape coefficients, in addition to complex amplitude in focal plane
        // the optimization code will then simultaneously minimize the light in the focal plane AND the PIAA shape deviations from the nominal shape
        if(REGPIAASHAPES==1)
        {
            // there are 4 groups of PIAA shape parameters: 2 sets of cosine modes and 2 sets of Fourier modes
            size1Dvec += data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1FmodesID].md[0].size[0];
        }

        // The same approach is used for regularization of the focal plane mask sag values
        // the sag values are appended to the evaluation vector
        if(REGFPMSAG==1)
        {
            size1Dvec += data.image[piaacmc[0].zonezID].md[0].size[0];
        }


        // re-package vector into 1D array and add regularization terms
        // the resulting vector is stored in vecDHref1D
        // we also create a mask image which is intended to mask out some of the pixels from the evaluation
        // the mask is currently not used, so we will write 1.0 in all of its pixels
        // DHmask and vecDHref1D contain all the optimization parameters as set above
        IDm = create_2Dimage_ID("DHmask", size1Dvec, 1); // "ID of evaluation mode mask"
        ID1Dref = create_2Dimage_ID("vecDHref1D", size1Dvec, 1); // "ID of 1D dark zone reference"

        // we first write 1.0 into the focal plane complex amplitudes in the vector
        ID = image_ID("vecDHref");
        for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
        {
            // imbed vecDHref into the evaluation zone part of of the full parameter vector vecDHref1D
            data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // sets the evaluation zone part of of the full parameter vector vecDHref1D to 1
            // 1 means the evaluation zone is on
            data.image[IDm].array.F[ii] = 1.0;
        }
        // !!!!!! WARNING !!!!!!!
        // the state of ii at this point drives the code below and will evolve until the comment
        // that says we're done with ii

        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 8;

        // Now actually fill in the regularized output vector.
        // If we are not regularizing, the output evaluation zone values are filled in by the
        // PSF calls in the optimization loop
        // and then append regularization vectors in the main evaluation vector
        // Note: index ii is incremented as we add "sub-vectors" into the main evaluation vector
        if(REGPIAASHAPES == 1)
        {
            // initialize by filling in the regularization terms of the output,
            // !! starting at the current value of ii !!
            // This means we're actually writing the output vector regularization terms.
            // Otherwise this is the same as the if(REGPIAASHAPES == 1) code block above

            ID = piaacmc[0].piaa0CmodesID;
            IDref = image_ID("piaa0Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID = piaacmc[0].piaa1CmodesID;
            IDref = image_ID("piaa1Cmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID_CPAfreq = image_ID("cpamodesfreq");
            ID = piaacmc[0].piaa0FmodesID;
            IDref = image_ID("piaa0Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }

            ID = piaacmc[0].piaa1FmodesID;
            IDref = image_ID("piaa1Fmref");
            if(IDref==-1)
            {
                IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                    data.image[IDref].array.F[jj] = 0.0;
            }
            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }




        // same for the sags, starting at the current value of ii
        if(REGFPMSAG == 1)
        {
            ID = piaacmc[0].zonezID;
            for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = pow( data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }
        // !!!!!!
        // we're done with ii


        delete_image_ID("vecDHref"); // vecDHref has beem embedded into vecDHref1D

        // at this point, we have completed the initialization, and the optimization loop starts


        initbestval = 0;
        // file that will track optimization loop progress
        sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "w");
        fclose(fp);

        //	list_image_ID();
        //
        // LINEAR OPTIMIZATION AROUND CURRENT POINT
        //

        iterOK=1;
        iter = 0;
        oldval = 1.0;
        data.image[IDstatus].array.UI16[0] = 9;

        // while # of iterations < NBiter
        //  and the ojective changes by more than 2% after the second iteration
        //  and something about NBlinoptgain ???????
        while(iterOK==1)//        for(iter=0; iter<NBiter; iter++)
        {
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 10;
            printf("Iteration %ld/%ld\n", iter, NBiter);
            fflush(stdout);
         
            // array for collecting dark hole mode derivatives
            // stores derivative of output vector against input parameters
            IDmodes = create_3Dimage_ID("DHmodes", size1Dvec, 1, number_param);
            // 2D array for diagnostic display
            IDmodes2D = create_2Dimage_ID("DHmodes2D", size1Dvec, number_param); //TEST

            // get ready to update optimization tracking file
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### PIAACMC_FPM_FASTDERIVATIVES = %d\n", piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES);
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 11;

            printf("Compute local derivatives of output vector against input focal plane mask zones (sags)\n");
            fflush(stdout);
            
            if( piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES == 1) // TO BE USED ONLY FOR FOCAL PLANE MASK OPTIMIZATION
            { // this only happens in mode 13
                // the fast derivative mode only works for focal plane mask optimization, for which derivatives against sag values can be comptuted by simple rotation of pre-computed vectors from mode 11
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 12;
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                //				ID = create_2Dimage_ID("DHmodes2Dtest", size1Dvec, number_param);
            
				printf("Computing %ld derivatives ", (long) data.image[piaacmc[0].zonezID].md[0].size[0]);
				fflush(stdout);
                for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++) // loop over mask zones
                {
					printf(" %ld", mz);
                    fflush(stdout);
                    // actually compute the derivative
                    // fpmresp_array is results from mode 11
                    // from mode 13 above:
                    //      fpmresp_array = data.image[IDfpmresp].array.D;
                    //      zonez_array = data.image[piaacmc[0].zonezID].array.D;
                    // dphadz_array was computed in mode 13 shortly afterwards
                    // outtmp_array is output
                    PIAACMCsimul_achromFPMsol_eval_zonezderivative(mz, piaacmcsimul_var.fpmresp_array, piaacmcsimul_var.zonez_array, piaacmcsimul_var.dphadz_array, piaacmcsimul_var.outtmp_array, piaacmcsimul_var.vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
                    for(ii=0; ii<size1Dvec0; ii++)
                        data.image[IDmodes].array.F[mz*size1Dvec+ii] = piaacmcsimul_var.outtmp_array[ii]*paramdelta[mz];
                }
                // derivatives of regularization values against sag values can also be computed analytically without requiring diffraction propagation
                if(REGFPMSAG == 1)
                {
                    ID = piaacmc[0].zonezID;
                    // following should be derivative of (sag/coeff)^alpha
                    // w.r.t. sag
                    for(mz=0; mz < data.image[ID].md[0].size[0]; mz++)
                        data.image[IDmodes].array.F[mz*size1Dvec + (size1Dvec0+mz)] = (piaacmc[0].fpmsagreg_alpha/piaacmc[0].fpmsagreg_coeff) * pow( data.image[ID].array.D[mz]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha-1.0)*paramdelta[mz];
                }


                // TEST diagnostic
                memcpy(data.image[IDmodes2D].array.F, data.image[IDmodes].array.F, sizeof(float)*size1Dvec*number_param);
                save_fl_fits("DHmodes2D", "!test_DHmodes2D.fits");

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 13;
            
				printf("Done computing derivatives (FAST MODE)\n");
				fflush(stdout);
            }
            else // ONLY FOR PIAA SHAPES OPTIMIZATION
            {
                // derivatives against PIAA shapes must be computed numerically
                
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 14;
                for(i=0; i<number_param; i++)
                {
                    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                    // get delta on-axis PSF response to the change in paramdelta
                    // to later give derivative w.r.t. paramdelta
                    if(paramtype[i]==_DATATYPE_FLOAT)
                    {
                        *(paramvalf[i]) += (float) paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    }
                    else
                    {
                        *(paramval[i]) += paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    }

                    //      sprintf(fname,"!%s/imvect_%02ld.fits", piaacmcsimul_var.piaacmcconfdir, i);
                    //       save_fits("imvect", fname);
                    ID = image_ID("imvect");


                    sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
                    fp = fopen(fname, "a");
                    fprintf(fp, "# %5ld/%5ld %5ld/%5ld %20.15g %20.15g %20.15g      %20.15g\n", iter, NBiter, i, number_param, paramdelta[i], val, valref, bestval);
                    fclose(fp);

                    // re-package vector into 1D array and add regularization terms
					// evaluation vector is "imvect1D", ID = ID1D
                    // similar to vecDHref before
                    ID1D = create_2Dimage_ID("imvect1D", size1Dvec, 1);
                    // fill in the evaluation point portion
                    for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                        data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii];

                    if(REGPIAASHAPES==1)
                    { // fill in the shape regularization value's response to paramdelta using the
                        // same formulas as before with the delta PSF as input
                        ID = piaacmc[0].piaa0CmodesID;
                        IDref = image_ID("piaa0Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa0C_regcoeff_alpha);
                            ii++;
                        }

                        ID = piaacmc[0].piaa1CmodesID;
                        IDref = image_ID("piaa1Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa1C_regcoeff_alpha);
                            ii++;
                        }

                        ID_CPAfreq = image_ID("cpamodesfreq");

                        ID = piaacmc[0].piaa0FmodesID;
                        IDref = image_ID("piaa0Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                            ii++;
                        }

                        ID = piaacmc[0].piaa1FmodesID;
                        IDref = image_ID("piaa1Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        {
                            data.image[ID1D].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                            ii++;
                        }
                    }


                    delete_image_ID("imvect"); // has been imbedded into imvect1D

					// restore original state (return to original staring point)
                    if(paramtype[i]==_DATATYPE_FLOAT)
                        *(paramvalf[i]) -= (float) paramdelta[i];
                    else
                        *(paramval[i]) -= paramdelta[i];


                    // compute actual derivative as first difference from reference
                    // this is the starting derivative
                    for(ii=0; ii<data.image[ID1D].md[0].nelement; ii++)
                        data.image[IDmodes].array.F[i*data.image[ID1D].md[0].nelement+ii] = (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);


                    //    printf("%3ld %g %g\n", i, val, valref);

                    // create diagnostic image
                    ID = create_2Dimage_ID("DHmodes2D", size1Dvec, number_param);
                    for(ii=0; ii<data.image[IDmodes].md[0].nelement; ii++)
                        data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];

                    sprintf(fname, "!%s/DMmodes.fits", piaacmcsimul_var.piaacmcconfdir);
                    save_fits("DHmodes2D", fname);

                    delete_image_ID("DHmodes2D");
                }
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 15;
            }


            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 16;
            // print the results to file for human tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### scanning gain \n");
            fprintf(fp, "### <alphareg>  <gain>  <contrast>\n");
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 17;

            // first three arguments are names of the input arrays
            // vecDHref1D is the input data
            // DHmodes is the basis of modes to expand vecDHref1D into
            // DHmask weights elements of vecDHref1D (nominally all weights = 1)
            // 4th arg is the pseudoinverse eigenvalue (via eigenvalue decomposition)
            // this decomposes vecDHref1D into the DHmodes
            // 5th arg is the output: optcoeff0*DHmodes = vecDHref1D
            // computed via pseudoinverse of DHmodes
            // This decomposition facilitates the cancellation of vecDHref1D by
            // searching in the DHmodes basis
            //
            // use three cutoff values to give three options for future evaluation
            // smallest cutoff values produce the largest changes (are least well conditioned)
            printf("ref = 0.1   -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.1, "optcoeff0", 0);
            printf("- DONE\n");
            fflush(stdout);            
            
            printf("ref = 0.01  -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.01, "optcoeff1", 0);
            printf("- DONE\n");
            fflush(stdout);    
            
            printf("ref = 0.001 -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.001, "optcoeff2", 0);
            printf("- DONE\n");
            fflush(stdout);                
            
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 18;

            // initialize zero "optimal" vector optvec giving direction to move in the search for the min
            arith_image_cstmult("optcoeff0", 0.0, "optvec"); // create optimal vector
            IDoptvec = image_ID("optvec");
            // initialize the objective value
            initbestval = 0;
            bestval = valref;

            alphareg = 1.0; // has no effect (see next loop)

            // say something here ???
            NBlinoptgain = 0;


            scangainfact = 1.2;
            alphascaninit = 0;
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 19;
            // alphareg controls linear combinations of the directions optcoeff0,1,2
            // alphareg = 0 => moving along optcoeff0
            // alphareg = 0.5 => moving along optcoeff1
            // alphareg = 1 => moving along optcoeff2
            // with interpolation
            // look in 5 steps if alphareg += 0.2
            for(alphareg=0.0; alphareg<1.01; alphareg += 0.2)
            {
                // produce piecewise linear interoplation coefficients
                acoeff0 = 1.0 - 2.0*alphareg; // ranges from -1 to 1
                if(acoeff0<0.0)
                    acoeff0 = 0.0; // clip to 0 to 1, = 0 if alphareg > 0.5

                acoeff1 = 1.0 - fabs(2.0*(alphareg-0.5)); // two lines: from 0,1 to .5,0 to 1,1

                acoeff2 = 2.0*alphareg - 1.0; // ranges from -1 to 1
                if(acoeff2<0.0)
                    acoeff2 = 0.0; // clip to 0 to 1, = 0 if alphareg < 0.5
                
                // sum of acoeff0,1,2 = 1 at all values of alphareg.

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 20;
                
                // optcoeff0m = acoeff0*optcoeff0 etc.
                arith_image_cstmult("optcoeff0", acoeff0, "optcoeff0m");
                arith_image_cstmult("optcoeff1", acoeff1, "optcoeff1m");
                arith_image_cstmult("optcoeff2", acoeff2, "optcoeff2m");

                // diagnostic
                save_fl_fits("optcoeff0", "!optcoeff0.fits");//TEST
                save_fl_fits("optcoeff1", "!optcoeff1.fits");
                save_fl_fits("optcoeff2", "!optcoeff2.fits");

                // optcoeff01m = acoeff0*optcoeff0 + acoeff1*optcoeff1
                arith_image_add("optcoeff0m", "optcoeff1m", "optcoeff01m");
                // optcoeff = acoeff0*optcoeff0 + acoeff1*optcoeff1 + acoeff2*optcoeff2
                arith_image_add("optcoeff01m", "optcoeff2m", "optcoeff");
                // optcoeff now has our search direction
                delete_image_ID("optcoeff0m");
                delete_image_ID("optcoeff1m");
                delete_image_ID("optcoeff2m");
                delete_image_ID("optcoeff01m");

                ID = image_ID("optcoeff");
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 21;

                // do linear scan along the direction optcoeff from current parameter location
                linscanOK = 1;
                // size of step in this direction
                scangain = 0.0; //scanstepgain; overriden below
                val = 100000000000.0; // initialize minimization objective to a big number
                bestgain = 0.0;
                // iteration counter of steps in the current direction optcoeff
                k = 0;

                // if(alphascaninit==1)
                //		scangain = bestgain/scangainfact/scangainfact/scangainfact/scangainfact/scangainfact;
                //	alphascaninit = 1;

                // scangain is our location along the direction optcoeff
                scangain = 0.001;
                //              scanstepgain = 0.000001; // TEST
                //              scangainfact = 1.00001; // TEST

                // while objective value < previous value and we've taken no more than than 90 steps
                while(linscanOK==1)
                {
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 22;

                    // compute offsets
                    ID = image_ID("optcoeff"); // direction vector
                    linoptlimflagarray[k] = 0;
                    // step each parameter by optcoeff
                    for(i=0; i<number_param; i++) // looping over parameters
                    {
                        // compute step delta for this parameter
                        // image[ID] = optcoeff is a derivative w.r.t. paramdelta, so has dimension
                        // (parameter dimension)/(paramdelta dimension) so we have to mulitply by
                        // paramdelta to put our step in physical parameter units
                        // negative because we want to cancel the value from the delta PSF
                        paramdeltaval[i] = -scangain * data.image[ID].array.F[i] * paramdelta[i];
                        if(paramdeltaval[i]<-parammaxstep[i]) // if the step is too large in the negative direction
                        {
                            printf("MIN LIMIT [%3ld   %20g]   %20g -> ", i, paramdelta[i], paramdeltaval[i]); //TEST
                            paramdeltaval[i] = -parammaxstep[i]; // set it to the negative largest allowed step
                            printf(" %20g\n", paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }
                        if(paramdeltaval[i]>parammaxstep[i])// if the step is too large in the positive direction
                        {
                            printf("MAX LIMIT [%3ld   %20g]   %20g -> ", i, paramdelta[i], paramdeltaval[i]); //TEST
                            paramdeltaval[i] = parammaxstep[i]; // set it to the positive largest allowed step
                            printf(" %20g\n", paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }

                        // apply offsets to the global data object via the pointers paramvalf, which
                        // point into the (hopefully) coorect locations of each parameter's value in the
                        // data object
                        if(paramtype[i]==_DATATYPE_FLOAT)
                        {
                            if(  *(paramvalf[i]) + (float) paramdeltaval[i]  > parammax[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                paramdeltaval[i] = parammax[i] - *(paramvalf[i]);

                            if(  *(paramvalf[i]) + (float) paramdeltaval[i]  < parammin[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                paramdeltaval[i] = parammin[i] - *(paramvalf[i]);
                            // take the actual step (paramvalf is a 1D array of pointers)
                            *(paramvalf[i]) += (float) paramdeltaval[i];
                        }
                        else // same for the double case
                        {
                            if(  *(paramval[i]) + paramdeltaval[i]  > parammax[i] )
                                paramdeltaval[i] = parammax[i] - *(paramval[i]);

                            if(  *(paramval[i]) + paramdeltaval[i]  < parammin[i] )
                                paramdeltaval[i] = parammin[i] - *(paramval[i]);

                            *(paramval[i]) += paramdeltaval[i];
                        }
                    }
                    // store the current objective value for later comparison
                    valold = val;
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 23;

					// compute new state and compute assossiated evaluation metric
                    // using the modified global data object
                    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    valContrast = val; // contrast component of the evaluation metric
                    // we've now only done the light portion
                    
                    // add regularization component of the evaluation metrix
                    // first, compute and add PIAA shape regularization value (val0) if applicable
                    // this is the same as the previous computation of val0 around line 7627
                    val0 = 1.0;
                    if(REGPIAASHAPES==1)
                    {
                        val0 = 0.0;
                        ID = piaacmc[0].piaa0CmodesID;
                        IDref = image_ID("piaa0Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }

                        ID = piaacmc[0].piaa1CmodesID;
                        IDref = image_ID("piaa1Cmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }


                        ID_CPAfreq = image_ID("cpamodesfreq");

                        ID = piaacmc[0].piaa0FmodesID;
                        IDref = image_ID("piaa0Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }

                        ID = piaacmc[0].piaa1FmodesID;
                        IDref = image_ID("piaa1Fmref");
                        if(IDref==-1)
                        {
                            IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                            for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                                data.image[IDref].array.F[jj] = 0.0;
                        }
                        for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        {
                            tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                            val0 += tmp*tmp;
                        }
                        val += val0;
                    }

					// add sag regularization (val1) if applicable, as before
                    val1 = 1.0;
                    if(REGFPMSAG == 1)
                    {
                        ID = piaacmc[0].zonezID;
                        val1 = 0.0;
                        for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                        {
                            tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                            val1 += tmp*tmp;
                        }
                        val += val1;
                    }
                    // val is now our complete objective!! Yay!!

                    
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 24;
                    // print it for monitoring
                    sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
                    fp = fopen(fname, "a");
                    // printf the first part of the line reporting current values
                    fprintf(fp, "##  %5.3f   %20lf           %20g    (reg = %12g %12g   contrast = %20g)       [%d] [%ld]", alphareg, scangain, val, val0, val1, valContrast, linoptlimflagarray[k], number_param);
 
                    // now add text indicating status and complete line
                    // and store all parameters for the current best solution
                    if((val<bestval)||(initbestval==0))
                    {
                        for(i=0; i<number_param; i++)
                            if(paramtype[i]==_DATATYPE_FLOAT)
                                data.image[IDoptvec].array.F[i] = *(paramvalf[i]);
                            else
                                data.image[IDoptvec].array.F[i] = (float) *(paramval[i]);
                        bestval = val;
                        if(initbestval == 0)
                            fprintf(fp, " ===== START POINT =====\n");
                        else
                            fprintf(fp, "  -> BEST VECTOR =======\n");
                        bestgain = scangain;
                        initbestval = 1;
                    }
                    else
                    {
                        fprintf(fp, " bestval = %12g\n", bestval);
                    }
                    fclose(fp);

                    // remove offsets returning the global data object to its original state
                    for(i=0; i<number_param; i++)
                    {
                        if(paramtype[i]==_DATATYPE_FLOAT)
                            *(paramvalf[i]) -= (float) paramdeltaval[i];
                        else
                            *(paramval[i]) -= paramdeltaval[i];
                    }
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 25;

                    // store the current position and value
                    linoptgainarray[k] = scangain;
                    linoptvalarray[k] = val;
                    k++; // next step

                    // test to see if we're no longer getting better
                    if(val<valold)
                    {
                        linscanOK = 1; // if we're getting better keep going
                        // bestgain = scangain;
                        // scangain += scanstepgain;
                    }
                    else // otherwise stop stepping
                        linscanOK = 0;



                    if(k>90)  // stop if we've taken too many steps
                        linscanOK = 0;

                    // increment our location along the line
                    scangain += scanstepgain; // scanstepgain is an initilizaed function local (currently 0.001)
                    scangain *= scangainfact; // causes later steps to be larger
                                                // (implicit scangainfact^n for the nth step)
                }
                // NBlinoptgain is counting the largest number of steps needed in this inner loop
                // stepping in the current direction.  When this is small (< 3) we declare victory
                // and stop the outer linear optimization
                if(k>NBlinoptgain)
                    NBlinoptgain = k;

                delete_image_ID("optcoeff"); // delete the current direction
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 26;
            }
            // best solution after this linear linescan is stored in IDoptvec
            delete_image_ID("optcoeff0");
            delete_image_ID("optcoeff1");
            delete_image_ID("DHmodes");

            // we've now found the minimum using the three directions from the
            // alternative decompositions of the parameter space with the DHmodes basis
            // (with different conditioning).
            // Now we check the result by recomputing from scratch the objective with the current
            // (hopefully) optimal parameters.  Belts and suspenders
            
            // At this point the global data object has been restored to its original
            // (non-optimal) state.
            
            // (re-)compute best solution identified in previous linescan
			// update state to best solution (IDoptvec), setting the global data object
            // to the optimal state
            for(i=0; i<number_param; i++)
            {
                if(paramtype[i]==_DATATYPE_FLOAT)
                    *(paramvalf[i]) = data.image[IDoptvec].array.F[i];
                else
                    *(paramval[i]) = (double) data.image[IDoptvec].array.F[i];
            }
            valold = val;
			// compute contrast metric component -> val using the data object in the latest optimal state
            val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);

			// add PIAA shape regularization component (val0) if applicable
            // same val0 and val1 computations are before, using the latest optimal state
            val0 = 1.0;
            if(REGPIAASHAPES==1)
            {
                val0 = 0.0;
                ID = piaacmc[0].piaa0CmodesID;
                IDref = image_ID("piaa0Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa0C_regcoeff_alpha);
                    val0 += tmp*tmp;
                }

                ID = piaacmc[0].piaa1CmodesID;
                IDref = image_ID("piaa1Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaa1C_regcoeff_alpha);
                    val0 += tmp*tmp;
                }


                ID_CPAfreq = image_ID("cpamodesfreq");

                ID = piaacmc[0].piaa0FmodesID;
                IDref = image_ID("piaa0Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                    val0 += tmp*tmp;
                }

                ID = piaacmc[0].piaa1FmodesID;
                IDref = image_ID("piaa1Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                {
                    tmp = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                    val0 += tmp*tmp;
                }
                val += val0;
            }
            
            // add sag regularization component if applicable
            val1 = 1.0;
            if(REGFPMSAG == 1)
            {
                ID = piaacmc[0].zonezID;
                val1 = 0.0;
                for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    tmp = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                    val1 += tmp*tmp;
                }
                val += val1;
            }
            
            // now val is the objective including any desired regularization terms using the
            // latest optimal solution


            printf("gain: %lf -> val = %20g\n", bestgain, val);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 27;



			// update reference state evaluation vector of optimal results including evaluation zone values
            // and desired regularization terms
            // which sets up starting the next iteration at the best solution
            
            ID1Dref = image_ID("vecDHref1D");
            ID = image_ID("imvect");
            // first fill in evaluation zone (complex) values
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 28;

            // now fill in the regularization terms if desired
            // same code as before.
            if(REGPIAASHAPES==1)
            {
                ID = piaacmc[0].piaa0CmodesID;
                IDref = image_ID("piaa0Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }

                for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa0C_regcoeff_alpha);
                    ii++;
                }

                ID = piaacmc[0].piaa1CmodesID;
                IDref = image_ID("piaa1Cmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj,piaa1C_regcoeff_alpha);
                    ii++;
                }


                ID_CPAfreq = image_ID("cpamodesfreq");

                ID = piaacmc[0].piaa0FmodesID;
                IDref = image_ID("piaa0Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) *pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa0F_regcoeff_alpha);
                    ii++;
                }

                ID = piaacmc[0].piaa1FmodesID;
                IDref = image_ID("piaa1Fmref");
                if(IDref==-1)
                {
                    IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
                    for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                        data.image[IDref].array.F[jj] = 0.0;
                }
                for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) *pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaa1F_regcoeff_alpha);
                    ii++;
                }
            }


            if(REGFPMSAG == 1)
            {
                ID = piaacmc[0].zonezID;
                for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                    ii++;
                }
            }

            delete_image_ID("imvect");
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 29;


            

       




            // print out current best value for tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            if(fp==NULL)
            {
                printf("ERROR: cannot open file \"%s\"\n", fname);
                exit(0);
            }
            fprintf(fp, "-> %5ld    %20g <- %20g \n", iter, val, valref);
            printf("%5ld %20g %20g \n", iter, val, valref);
            fflush(stdout);
            fclose(fp);

            // save current best value and reference value in globals
            piaacmcsimul_var.PIAACMCSIMUL_VAL = val;
            piaacmcsimul_var.PIAACMCSIMUL_VALREF = valref;



            // Nominally if we're in this linear
            // optimization piaacmcsimul_var.PIAACMC_fpmtype = 1, so the next line is not executed
            if(piaacmcsimul_var.PIAACMC_fpmtype==0) // in the idealized PIAACMC case
                piaacmc[0].fpmaskamptransm = data.image[piaacmc[0].zoneaID].array.D[0]; // required to ensure that the new optimal focal plane mask transmission is written to disk

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 30;


            // tracking diagnostics giving behavior of the modes by iteration
            sprintf(dirname, "%s_linopt", piaacmcsimul_var.piaacmcconfdir);
            PIAACMCsimul_savepiaacmcconf(dirname); // staging area
            sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcsimul_var.piaacmcconfdir);
            r = system(command);

            sprintf(command, "cp %s/piaa0Cmodes.fits %s/piaa0Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa0Fmodes.fits %s/piaa0Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaa1Cmodes.fits %s/piaa1Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa1Fmodes.fits %s/piaa1Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaacmcparams.conf %s/piaacmcparams.%04ld.conf", dirname, dirname, iter);
            r = system(command);

			if(file_exists(stopfile)==1)
				iterOK = 0;

            // Figure out if current loop should continue optimization
            // if optimization ends, then iterOK set to 0
            // if we've reached the allowed number of iterations
            if(iter==NBiter)
                iterOK = 0;
            if(iter>2)
            {
                if(val>0.98*oldval) // if after second iteration and we'be improved by less than 10%
                    iterOK = 0;
            }
            
            if(NBlinoptgain<3) // if we've stopped moving much
                iterOK = 0;


            // set up for next iteration
            oldval = val;
            iter++;

            printf("END OF LOOP ITERATION\n");
            fflush(stdout);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 31;
        }
        printf(" ============ END OF OPTIMIZATION LOOP ======= \n");
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 32;
    } // end of if (piaacmcsimul_var.LINOPT==1): done with the linear optimization

	piaacmcsimul_var.LINOPT = 0;



    //  PIAACMCsimul_savepiaacmcconf("piaacmc0");
    //  PIAACMCsimul_loadpiaacmcconf("piaacmc0");
    // PIAACMCsimul_savepiaacmcconf("piaacmc1");
    //exit(0);

    return 0;
}


