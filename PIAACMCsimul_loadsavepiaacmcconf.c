/**
 * @file    PIAAACMCsimul_loadsavepiaacmcconf.c
 * @brief   PIAA-type coronagraph design
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



// milk includes
//   core modules
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
//   other modules



#include "PIAACMCsimul/PIAACMCsimul.h"





extern DATA data;   

extern PIAACMCsimul_var piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;




int PIAAsimul_loadpiaacmcconf(const char *dname)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[500];
    char imname[500];
    long i;

    int tmpi;
    long tmpl;
    float tmpf;
    double tmplf;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    sprintf(fname,"%s/piaacmcparams.conf", dname);
    printf("%s\n", fname);


    fp = fopen(fname, "r");
    if(fp==NULL)
    {
        printf("Configuration file \"%s\" does not exist (yet), using previously set configuration\n", fname);
        fflush(stdout);
        r = 0;
    }
    else
    {

        r = fscanf(fp, "%ld   NBradpts\n", &tmpl);
        piaacmc[0].NBradpts = tmpl;

         for(i=0; i<10; i++)
        {
            if(i<piaacmc[0].NBLyotStop)
            {
                sprintf(fname, "%s/LyotStop%ld.fits", dname, i);
                sprintf(imname, "lyotstop%ld", i);
                printf("Loading \"%s\" as \"%s\"\n", fname, imname);
                piaacmc[0].IDLyotStop[i] = load_fits(fname, imname, 1);
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
            else
            {
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
        printf("LYOT STOP %ld POS : %lf\n", i, tmplf);
            
        }


        sprintf(fname, "%s/piaa0Cmodes.fits", dname);
        piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        if(piaacmc[0].piaa0CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", dname);
            piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        }


        sprintf(fname, "%s/piaa0Fmodes.fits", dname);
        piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        if(piaacmc[0].piaa0FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", dname);
            piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Cmodes.fits", dname);
        piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        if(piaacmc[0].piaa1CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", dname);
            piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Fmodes.fits", dname);
        piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        if(piaacmc[0].piaa1FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", dname);
            piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        }


        r = fscanf(fp, "%f   fpmaskamptransm\n",    &tmpf);
        piaacmc[0].fpmaskamptransm = tmpf;

        r = 1;

        fclose(fp);
    }
    
    return(r);
}







int PIAAsimul_savepiaacmcconf(const char *dname)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[500];
    long i;


	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif

    
    sprintf(command, "mkdir -p %s", dname);
    r = system(command);

    sprintf(fname,"%s/piaacmcparams.conf", dname);
    fp = fopen(fname, "w");


    fprintf(fp, "%10ld   NBradpts\n", piaacmc[0].NBradpts);

    for(i=0; i<10; i++)
    {
        if(i<piaacmc[0].NBLyotStop)
        {
            sprintf(fname, "!%s/LyotStop%ld.fits", dname, i);
            if(piaacmc[0].IDLyotStop[i]!=-1)
                save_fits(data.image[piaacmc[0].IDLyotStop[i]].name, fname);
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
        }
        else
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
    }

    sprintf(fname, "!%s/piaa0Cmodes.fits", dname);
    if(piaacmc[0].piaa0CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0CmodesID].name, fname);

    sprintf(fname, "!%s/piaa0Fmodes.fits", dname);
    if(piaacmc[0].piaa0FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0FmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Cmodes.fits", dname);
    if(piaacmc[0].piaa1CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1CmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Fmodes.fits", dname);
    if(piaacmc[0].piaa1FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1FmodesID].name, fname);


    sprintf(fname, "!%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

    
    if(piaacmc[0].zonezID!=-1)
        save_fits(data.image[piaacmc[0].zonezID].name, fname);

    fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmc[0].fpmaskamptransm);

    sprintf(fname, "!%s/fpm_zonea_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_minsag%06ld_maxsag%06ld_fpmreg%06ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1000.0*piaacmc[0].fpmsagreg_coeff+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

    
    if(piaacmc[0].zoneaID!=-1)
        save_fits(data.image[piaacmc[0].zoneaID].name, fname);

    fclose(fp);

    

    return(0);
}


