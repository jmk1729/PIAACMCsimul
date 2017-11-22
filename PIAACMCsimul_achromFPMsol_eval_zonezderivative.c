/**
 * @file    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
 * @brief   PIAA-type coronagraph design, run
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



#include <stdlib.h>
#include <stdio.h>
#include <math.h>












double PIAACMCsimul_achromFPMsol_eval_zonezderivative(long zone, double *fpmresp_array, double *zonez_array, double *dphadz_array, double *outtmp_array, long vsize, long nbz, long nbl)
{
	long evalk;

    // axis 0: eval pts (ii)   size = data.image[IDfpmresp].md[0].size[0] -> vsize
    // axis 1: zones (mz)      size = data.image[piaacmc[0].zonezID].md[0].size[0]+1 = nbz+1
    // axis 3: lambda (k)      size = piaacmc[0].nblambda -> nbl
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii

	#ifdef PIAASIMUL_LOGFUNC1
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif

	long evalki1, evalmz, evalki, evalkv;
	double evalpha, evalcosp, evalsinp;
    for(evalk=0; evalk<nbl; evalk++) // lambda loop
    {
        evalki = evalk*(nbz+1)*vsize;


        // outer zone
        //  for(evalii=0; evalii<vsize; evalii++)
        //    outtmp_array[evalk*vsize+evalii] = 0.0; //fpmresp_array[evalk*(nbz+1)*vsize+evalii];


        evalmz = zone;
        
        // compute derivative as
        // dphadz is a function of wavelength
        // zonez is the current zone thickness
        // -zonez*dphadz sets the phase gives the resulting phase
        // Re = Re * -dphadz*sin(-zonez*dphadz) - Im * dphadz*cos(-zonez*dphadz)
        // Im = Re * dphadz*cos(-zonez*dphadz) + Im * -dphadz*sin(-zonez*dphadz)

        evalpha = -zonez_array[evalmz]*dphadz_array[evalk];
        // !!! note that cos is sin and sin is cos !!!
        // this implements a 90 degree pre-rotation so that this is a
        // derivative
        evalcosp = sin(evalpha)*dphadz_array[evalk]; //cos(evalpha);
        evalsinp = -cos(evalpha)*dphadz_array[evalk]; //sin(evalpha);
        evalki1 = evalki + (evalmz+1)*vsize;
        evalkv = evalk*vsize;

		long evalii;
		long evalii1, evalii2;
		double evalre, evalim, evalre1, evalim1;
        for(evalii=0; evalii<vsize/2; evalii++)
        {
            evalii1 = 2*evalii;
            evalii2 = 2*evalii+1;
            evalre = fpmresp_array[evalki1 + evalii1];
            evalim = fpmresp_array[evalki1 + evalii2];
            evalre1 = evalre*evalcosp - evalim*evalsinp;
            evalim1 = evalre*evalsinp + evalim*evalcosp;
            outtmp_array[evalkv + evalii1] = evalre1;
            outtmp_array[evalkv + evalii2] = evalim1;
        }
    }

    return 0.0;
}


