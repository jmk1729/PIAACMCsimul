
# Design steps: Monochromatic design ( #step < 100 ) {#page_PIAACMCsimul_DesignSteps_Mono}




## Overview


The design proceeds in discrete steps. The example scripts show the individual steps, and each step is executed with a separate command line. 

Steps from 0 to 99 are executed sequentially to design a monochromatic PIAACMC. For these steps, the user may run multiple steps with a single command. For example, running step 18 will execute all steps from 0 to 17 included. If a step has already been completed, it will not be re-run.







The monochromatic PIAACMC design process is as follows:

* design an idealized monochromatic PIAACMC for a centrally obscured aperture (steps 1-4)

* modify the design for the pupil aperture (steps 5-)





## STEP 000 (MODE=0): Create an idealized centrally obscured apodized PIAACMC monochromatic design

This is meant as a starting point for the PIAACMC, which will then be optimized further. This step takes a few minutes, and upon normal completion, display: 

~~~
./runPIAACMC REACHED STATE EXIT POINT (1)
~~~


The prolate function for a centrally obscured pupili is created, along with the idealized focal plane mask. The diffraction propagation code then computes complex amplitude in each plane.

This step takes a few minutes - most of the time is spent on iterations to compute the 2D apodization prolate function, executing discrete Fourier transforms (DFTs).

After this step, the contrast will likely be around 1e-7 to 1e-8. The next step to improve this nominal design is to find the optimal locations for the Lyot stops.


Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
ref/Cmodes_1024.fits                    | Cosine modes applied for PIAA optics shapes optimization
ref/Fmodes_1024.fits                    | Fourier modes applied for PIAA optics shapes optimization
piaacmcconf_i000/APLCapo.<X>.<Y>.info   | Apodization function info
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/conjugations.txt       | List of planes and conjugation distance to reference
piaacmcconf_i000/apo2Drad.fits          | Amplitude apodization (entirely allocated to PIAA optics)
piaacmcconf_i000/PIAA_Mshapes.txt       | Aspheric optics shapes (r0 z0 r1 z1), unit [m]
piaacmcconf_i000/piaam0z.fits           | First PIAA mirror sag [m]
piaacmcconf_i000/piaam1z.fits           | Second PIAA mirror sag [m]
piaacmcconf_i000/psfi0_step000.fits     | Coronagraphic PSF (flux normalized to total=1 without coronagraph)
log/varlog.txt                          | List of coronagraph parameters



The conjugations.txt file should contain number, location and description of each plane:

~~~
00  0.000000    input pupil
01  0.000000    TT mirror
02  -1.102609   pupil plane apodizer
03  1.199997    PIAA optics 0
04  -1.102609   PIAA optics 1
05  -1.102609   opaque mask at PIAA elem 1
06  -1.102609   post focal plane mask pupil
07  0.000000    Lyot mask 0
08  0.000000    back end pupil stop  (rad = 0.920000)
~~~

The complex amplitude for each of these planes (files WFamp0_xxx.fits) should show the pupil apodization and coronagraphic effect induced by the focal plane mask which moves light into the central obstruction.



[Amplitude in each of the 9 planes](../src/PIAACMCsimul/doc/figures/examplePIAA_step0_WFamp.png)


The script provides some approximate performance metrics upon normal exit:

~~~
TOTAL = 0.002426
    FLUX   0    114332.0000 1.000000
    FLUX   1    114332.0000 1.000000
    FLUX   2    114331.9872 1.000000
    FLUX   3    114331.9668 1.000000
    FLUX   4    114331.9357 0.999999
    FLUX   5    114331.8634 0.999999
    FLUX   6     23794.1440 0.208114
    FLUX   7       551.2827 0.004822
    FLUX   8       277.3360 0.002426
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 26231.4 -> 2.00672e-06
optsyst[0].flux[0]  = 114332
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 8.58033e+06, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.69458e-07
~~~


Runtime 00:01:35


---


## STEP 001: Propagate solution and compute PSF

The previous solution is propagated. This should yield the same result as step 0. 

Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/psfi0_step001.fits     | On-axis PSF 


~~~
TOTAL = 0.002387
    FLUX   0    114332.0000 1.000000
    FLUX   1    114332.0000 1.000000
    FLUX   2    114331.9872 1.000000
    FLUX   3    114331.9668 1.000000
    FLUX   4    114331.9356 0.999999
    FLUX   5    114331.9105 0.999999
    FLUX   6     23785.0351 0.208035
    FLUX   7       544.8861 0.004766
    FLUX   8       272.9253 0.002387
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 25740.7 -> 1.96918e-06
optsyst[0].flux[0]  = 114332
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 8.18737e+06, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.57118e-07
saving contrast value [2.57118e-07] -> piaacmcconf_i000/contrast_ptsrc_sm0_s0_l0565_sr10_nbr001_mr000_minsag-10000_maxsag010000_fpmreg001000_ssr00_ssm0_Mirror_wb01.txt
~~~


Runtime 00:00:07

---

## STEP 002: Specify input pupil geometry

The pupil geometry is copied to file `piaacmcconf_i000/pupa0_1024.fits`

Runtime 00:00:00

---





## STEP 003 (mode = 0): compute on-axis PSF for new pupil geometry

Same as step 1, but taking into account the pupil geometry. With spiders or gaps in the pupil, the contrast will not be as good.

Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/psfi0_step003.fits     | On-axis PSF 


~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7614 1.000000
    FLUX   4    416183.7012 0.999999
    FLUX   5    155983.4386 0.374794
    FLUX   6     32825.0378 0.078871
    FLUX   7     30509.7138 0.073308
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 2.74249e+07 -> 0.000158334
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 7.35517e+09, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 1.74319e-05
~~~

Runtime 00:00:07

---


## STEP 004 (mode = 5): Compute Lyot stops shapes and locations, 1st pass

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search.


When the optimization completes, the best solutions are listed:

~~~
BEST SOLUTION: 0.000000000000 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.706241153417  0.019795686221
BEST SOLUTION: 6.695657648611 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.706237596354  0.019795151562
BEST SOLUTION: 7.365223413472 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.706219734865  0.019792653140
BEST SOLUTION: 8.034789178333 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.706026161578  0.019768034052
BEST SOLUTION: 8.704354943194 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.705220276320  0.019673042636
BEST SOLUTION: 9.373920708056 / 33.478288243056    0.000000556868 / 0.000044544711  -> 0.702039509211  0.019321456405
BEST SOLUTION: 32.808722478194 / 33.478288243056    0.000000668241 / 0.000044544711  -> 0.733651721375  0.019318291526
BEST SOLUTION: 33.478288243056 / 33.478288243056    0.000000668241 / 0.000044544711  -> 0.710514365016  0.018620791690
~~~

The optimal Lyot stop(s) conjugation(s) is written in file `piaacmcparams_step004.conf`

This step takes about 2hr.



## STEP 005 (mode = 2): Optimize focal plane mask transmission, 1st pass

See source code in PIAACMCsimul_exec()


Optimizes the focal plane mask transmission. Results are written in file `result_fpmt.log`

~~~
 0.031579 1.35944e-05  0 0.3 0.1
 0.131579 7.96891e-06  0 0.3 0.1
 0.231579 3.83775e-06  0 0.3 0.1
 0.331579 1.20092e-06  0 0.3 0.1
 0.431579 5.84203e-08  0 0.3 0.1
 0.531579 4.10247e-07  0 0.3 0.1
 0.631579 2.2564e-06  0 0.3 0.1

 0.341579 1.01943e-06  1 0.09 0.03
 0.371579 5.64602e-07  1 0.09 0.03
 0.401579 2.44266e-07  1 0.09 0.03
 0.431579 5.84203e-08  1 0.09 0.03
 0.461579 7.06368e-09  1 0.09 0.03
 0.491579 9.01969e-08  1 0.09 0.03
 0.521579 3.0782e-07  1 0.09 0.03

 0.434579 4.72326e-08  2 0.027 0.009
 0.443579 2.1739e-08  2 0.027 0.009
 0.452579 8.34926e-09  2 0.027 0.009
 0.461579 7.06368e-09  2 0.027 0.009
 0.470579 1.78822e-08  2 0.027 0.009
 0.479579 4.08046e-08  2 0.027 0.009
 0.488579 7.58315e-08  2 0.027 0.009

 0.453479 7.676e-09  3 0.0081 0.0027
 0.456179 6.38254e-09  3 0.0081 0.0027
 0.458879 6.17844e-09  3 0.0081 0.0027
 0.461579 7.06368e-09  3 0.0081 0.0027
 0.464279 9.03836e-09  3 0.0081 0.0027
 0.466979 1.21023e-08  3 0.0081 0.0027
 0.469679 1.62556e-08  3 0.0081 0.0027

 0.456449 6.31312e-09  4 0.00243 0.00081
 0.457259 6.17018e-09  4 0.00243 0.00081
 0.458069 6.12529e-09  4 0.00243 0.00081
 0.458879 6.17844e-09  4 0.00243 0.00081
 0.459689 6.32964e-09  4 0.00243 0.00081
 0.460499 6.57888e-09  4 0.00243 0.00081
 0.461309 6.92612e-09  4 0.00243 0.00081

 0.457340 6.16129e-09  5 0.000729 0.000243
 0.457583 6.14047e-09  5 0.000729 0.000243
 0.457826 6.12846e-09  5 0.000729 0.000243
 0.458069 6.12529e-09  5 0.000729 0.000243
 0.458312 6.13094e-09  5 0.000729 0.000243
 0.458555 6.14542e-09  5 0.000729 0.000243
 0.458798 6.16872e-09  5 0.000729 0.000243
~~~

This takes about 30mn


## STEP 006 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput

Takes about 10mn


## STEP 007 (mode = 40): Tune PIAA shapes and focal plane mask transm, 10 cosine modes, 5 Fourier modes

This takes approximately 20mn for size = 1024.

Progress can be tracked by watching file :

~~~
tail -f linoptval.txt
~~~


## STEP 008 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes


~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7608 1.000000
    FLUX   4    416183.6991 0.999999
    FLUX   5    189642.2892 0.455669
    FLUX   6       216.7150 0.000521
    FLUX   7       216.7150 0.000521
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 6272.06 -> 3.62109e-08
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 291196, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 6.90139e-10
~~~

This step takes 3 hr



## STEP 009 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput


This step takes approximately 4mn.



## STEP 010 (mode = 1): Tune Lyot stops conjugations

To view result:

~~~
tail -f result_LMpos.log
~~~

~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7608 1.000000
    FLUX   4    416183.6991 0.999999
    FLUX   5    189690.4898 0.455785
    FLUX   6       159.9713 0.000384
    FLUX   7       159.9713 0.000384
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 7922.02 -> 4.57368e-08
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 298404, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 7.07224e-10
~~~

This step takes approximately 5mn.


## STEP 011 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes

~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7605 1.000000
    FLUX   4    416183.6985 0.999999
    FLUX   5    192668.7345 0.462941
    FLUX   6       159.3355 0.000383
    FLUX   7       159.3355 0.000383
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 8636.87 -> 4.98638e-08
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 209897, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 4.97461e-10
~~~


## STEP 012 (mode = 40): Tune PIAA shapes and focal plane mask transm,  40 cosine modes, 150 Fourier modes

~~~
The total number of free parameters is 380 = (40+150)*2, so this routine takes a long time to complete (hours).
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7615 1.000000
    FLUX   4    416183.6991 0.999999
    FLUX   5    193239.3033 0.464312
    FLUX   6       159.4789 0.000383
    FLUX   7       159.4789 0.000383
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 306.006 -> 1.76669e-09
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 29496.3, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 6.99069e-11
~~~

takes 4hr

## STEP 013 (mode = 5): Compute Lyot stops shapes and locations, 3nd pass, 70% throughput

takes 4mn


## STEP 014 (mode = 1): Tune Lyot stops conjugations

~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7615 1.000000
    FLUX   4    416183.6991 0.999999
    FLUX   5    193287.4013 0.464428
    FLUX   6       176.5218 0.000424
    FLUX   7       176.5218 0.000424
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 2220.58 -> 1.28202e-08
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 147416, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 3.49379e-10
~~~

Takes 12mn

## STEP 015 (mode = 40): tune PIAA shapes and focal plane mask transm, 20 cosine modes, 20 Fourier modes

~~~
    FLUX   0    416183.9306 1.000000
    FLUX   1    416183.9306 1.000000
    FLUX   2    416183.8669 1.000000
    FLUX   3    416183.7606 1.000000
    FLUX   4    416183.6984 0.999999
    FLUX   5    194898.6729 0.468299
    FLUX   6       176.7694 0.000425
    FLUX   7       176.7694 0.000425
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
Peak constrast (rough estimate)= 1145.55 -> 6.6137e-09
optsyst[0].flux[0]  = 416184
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 86266.2, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.04453e-10
~~~

Takes 1.5hr

## STEP 016 (mode = 40): tune PIAA shapes and focal plane mask transm, 40 cosine modes, 150 Fourier modes
In this example step 17 is skipped by typing:

~~~
touch piaacmcconf_i000/step016.txt
~~~



## STEP 017 (mode = 40): tune PIAA shapes and focal plane mask transm, 40 cosine modes, 625 Fourier modes
In this example step 17 is skipped by typing:

~~~
touch piaacmcconf_i000/step017.txt
~~~
