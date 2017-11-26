
# Design steps: Monochromatic design ( #step < 100 ) {#page_PIAACMCsimul_DesignSteps_Mono}




## Overview


The design proceeds in discrete steps. The example scripts show the individual steps, and each step is executed with a separate command line. 

Steps from 0 to 99 are executed sequentially to design a monochromatic PIAACMC. For these steps, the user may run multiple steps with a single command. For example, running step 18 will execute all steps from 0 to 17 included. If a step has already been completed, it will not be re-run.







The monochromatic PIAACMC design process is as follows:

* design an idealized monochromatic PIAACMC for a centrally obscured aperture (steps 1-4)

* modify the design for the pupil aperture (steps 5-)


---

---



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


Runtime 00:01:36 (size 1024)


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

Runtime 00:00:00 (size 1024)

---





## STEP 003 (mode = 0): compute on-axis PSF for new pupil geometry

Same as step 1, but taking into account the pupil geometry. With spiders or gaps in the pupil, the contrast will not be as good.

Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/psfi0_step003.fits     | On-axis PSF 


~~~
TOTAL = 0.140551
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9359 0.999999
    FLUX   5    103003.6175 0.999996
    FLUX   6     27580.8543 0.267765
    FLUX   7     15763.2885 0.153036
    FLUX   8     14477.3372 0.140551
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 4.65338e+07 -> 0.00438592
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 3.77881e+09, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 0.000146207
~~~

Runtime 00:00:07 (size 1024)

---


## STEP 004 (mode = 5): Compute Lyot stops shapes and locations, 1st pass

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search.


When the optimization completes, the best solutions are listed:

~~~
BEST SOLUTION: 5.646381230377 / 12.274741805168    0.000000942068 / 0.000021030879  -> 0.700625751251  0.039432599420
BEST SOLUTION: 5.891876066481 / 12.274741805168    0.000000942068 / 0.000021030879  -> 0.700355617330  0.039385771373
BEST SOLUTION: 8.837814099721 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.721432709799  0.039346632150
BEST SOLUTION: 9.083308935825 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.716777237262  0.038827027960
BEST SOLUTION: 9.328803771928 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.712835211313  0.038398704677
BEST SOLUTION: 9.574298608031 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.708858067020  0.037977986529
BEST SOLUTION: 9.819793444135 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.705292934047  0.037609975658
BEST SOLUTION: 10.065288280238 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.701842146347  0.037262905586
BEST SOLUTION: 11.292762460755 / 12.274741805168    0.000001356577 / 0.000021030879  -> 0.703163640577  0.036896957877
BEST SOLUTION: 12.029246969065 / 12.274741805168    0.000001627893 / 0.000021030879  -> 0.701096376349  0.036559137596
~~~


The optimal Lyot stop(s) conjugation(s) is written in file `piaacmcconf_i000/piaacmcparams_step004.conf`

Runtime 00:53:49 (size 1024)

---


## STEP 005 (mode = 2): Optimize focal plane mask transmission, 1st pass

See source code in PIAACMCsimul_exec_optimize_fpmtransmission()


Optimizes the focal plane mask transmission. Results are written in file `piaacmcconf_i000/result_fpmt.log`

~~~
 -0.05940301  3.24056e-05         0          0.3          0.1
 +0.04059699  1.92755e-05         0          0.3          0.1
 +0.14059699  9.55978e-06         0          0.3          0.1
 +0.24059699  3.25834e-06         0          0.3          0.1
 +0.34059699  3.71207e-07         0          0.3          0.1
 +0.44059699  8.98387e-07         0          0.3          0.1
 +0.54059699  4.83988e-06         0          0.3          0.1

 +0.25059699  2.81598e-06         1         0.09         0.03
 +0.28059699  1.69377e-06         1         0.09         0.03
 +0.31059699  8.78845e-07         1         0.09         0.03
 +0.34059699  3.71207e-07         1         0.09         0.03
 +0.37059699  1.70859e-07         1         0.09         0.03
 +0.40059699  2.77799e-07         1         0.09         0.03
 +0.43059699  6.92026e-07         1         0.09         0.03

 +0.34359699  3.37345e-07         2        0.027        0.009
 +0.35259699  2.54193e-07         2        0.027        0.009
 +0.36159699  1.98699e-07         2        0.027        0.009
 +0.37059699  1.70859e-07         2        0.027        0.009
 +0.37959699  1.70676e-07         2        0.027        0.009
 +0.38859699  1.98148e-07         2        0.027        0.009
 +0.39759699  2.53276e-07         2        0.027        0.009

 +0.37149699  1.69596e-07         3       0.0081       0.0027
 +0.37419699  1.67467e-07         3       0.0081       0.0027
 +0.37689699  1.67827e-07         3       0.0081       0.0027
 +0.37959699  1.70676e-07         3       0.0081       0.0027
 +0.38229699  1.76014e-07         3       0.0081       0.0027
 +0.38499699  1.83841e-07         3       0.0081       0.0027
 +0.38769699  1.94156e-07         3       0.0081       0.0027

 +0.37176699  1.69271e-07         4      0.00243      0.00081
 +0.37257699  1.68446e-07         4      0.00243      0.00081
 +0.37338699  1.67844e-07         4      0.00243      0.00081
 +0.37419699  1.67467e-07         4      0.00243      0.00081
 +0.37500699  1.67314e-07         4      0.00243      0.00081
 +0.37581699  1.67384e-07         4      0.00243      0.00081
 +0.37662699  1.67679e-07         4      0.00243      0.00081

 +0.37427799  1.67442e-07         5     0.000729     0.000243
 +0.37452099  1.67379e-07         5     0.000729     0.000243
 +0.37476399  1.67336e-07         5     0.000729     0.000243
 +0.37500699  1.67314e-07         5     0.000729     0.000243
 +0.37524999  1.67311e-07         5     0.000729     0.000243
 +0.37549299  1.67329e-07         5     0.000729     0.000243
 +0.37573599  1.67367e-07         5     0.000729     0.000243
~~~

Runtime 00:03:56 (size 1024)


---


## STEP 006 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput

Source code: 

~~~
BEST SOLUTION: 30.006465758591 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.602147984170  0.001772430096
BEST SOLUTION: 31.006681283877 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.601618355971  0.001755040537
BEST SOLUTION: 32.006896809163 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.601141999997  0.001739932914
BEST SOLUTION: 33.007112334450 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.600760087273  0.001728177742
BEST SOLUTION: 34.007327859736 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.600397821167  0.001717338698
BEST SOLUTION: 44.009483112599 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.602141259615  0.001712595775
BEST SOLUTION: 45.009698637886 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.601238947002  0.001692309270
BEST SOLUTION: 46.009914163172 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.600243611414  0.001670430013
~~~



Runtime 00:02:50 (size 1024)




---


## STEP 007 (mode = 40): Tune PIAA shapes and focal plane mask transm, 10 cosine modes, 5 Fourier modes



Progress can be tracked by watching file :

~~~
tail -f linoptval.txt
~~~

The file shows contrast value improving. Lines starting with '##' are evaluations along lines of steepest gradient, and lines starting with `###` are evaluations of local derivatives.

Last lines of file `linoptval.txt`:

~~~
##  1.000               0.123419                    2.12854e-08    (reg =            1            1   contrast =          2.12854e-08)       [0] [31] bestval =  2.12103e-08
##  1.000               0.149303                    2.12684e-08    (reg =            1            1   contrast =          2.12684e-08)       [0] [31] bestval =  2.12103e-08
##  1.000               0.180363                    2.12484e-08    (reg =            1            1   contrast =          2.12484e-08)       [0] [31] bestval =  2.12103e-08
##  1.000               0.217636                    2.12264e-08    (reg =            1            1   contrast =          2.12264e-08)       [0] [31] bestval =  2.12103e-08
##  1.000               0.262363                    2.12032e-08    (reg =            1            1   contrast =          2.12032e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.316036                    2.11772e-08    (reg =            1            1   contrast =          2.11772e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.380443                    2.11525e-08    (reg =            1            1   contrast =          2.11525e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.457732                    2.11278e-08    (reg =            1            1   contrast =          2.11278e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.550478                    2.11099e-08    (reg =            1            1   contrast =          2.11099e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.661774                    2.11015e-08    (reg =            1            1   contrast =          2.11015e-08)       [0] [31]  -> BEST VECTOR =======
##  1.000               0.795329                     2.1113e-08    (reg =            1            1   contrast =           2.1113e-08)       [0] [31] bestval =  2.11015e-08
->     9             2.11015e-08 <-          2.14234e-07 
~~~

Runtime 01:30:31 (size 1024)


----


## STEP 008 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes


Progress can be tracked by watching file :

~~~
tail -f linoptval.txt
~~~

Last lines of file `linoptval.txt`:

~~~
### scanning gain 
### <alphareg>  <gain>  <contrast>
##  0.000               0.001000                    6.45963e-09    (reg =            1            1   contrast =          6.45963e-09)       [0] [81] ===== START POINT =====
##  0.000               0.002400                    6.45963e-09    (reg =            1            1   contrast =          6.45963e-09)       [0] [81] bestval =  6.45963e-09
##  0.200               0.001000                    6.45965e-09    (reg =            1            1   contrast =          6.45965e-09)       [0] [81] bestval =  6.45963e-09
##  0.200               0.002400                    6.45968e-09    (reg =            1            1   contrast =          6.45968e-09)       [0] [81] bestval =  6.45963e-09
##  0.400               0.001000                    6.45972e-09    (reg =            1            1   contrast =          6.45972e-09)       [0] [81] bestval =  6.45963e-09
##  0.400               0.002400                    6.45991e-09    (reg =            1            1   contrast =          6.45991e-09)       [0] [81] bestval =  6.45963e-09
##  0.600               0.001000                    6.45944e-09    (reg =            1            1   contrast =          6.45944e-09)       [0] [81]  -> BEST VECTOR =======
##  0.600               0.002400                    6.45992e-09    (reg =            1            1   contrast =          6.45992e-09)       [0] [81] bestval =  6.45944e-09
##  0.800               0.001000                    6.45954e-09    (reg =            1            1   contrast =          6.45954e-09)       [0] [81] bestval =  6.45944e-09
##  0.800               0.002400                    6.45906e-09    (reg =            1            1   contrast =          6.45906e-09)       [0] [81]  -> BEST VECTOR =======
##  0.800               0.004080                    6.45901e-09    (reg =            1            1   contrast =          6.45901e-09)       [0] [81]  -> BEST VECTOR =======
##  0.800               0.006096                    6.45865e-09    (reg =            1            1   contrast =          6.45865e-09)       [0] [81]  -> BEST VECTOR =======
##  0.800               0.008515                    6.45912e-09    (reg =            1            1   contrast =          6.45912e-09)       [0] [81] bestval =  6.45865e-09
##  1.000               0.001000                    6.45904e-09    (reg =            1            1   contrast =          6.45904e-09)       [0] [81] bestval =  6.45865e-09
##  1.000               0.002400                    6.45889e-09    (reg =            1            1   contrast =          6.45889e-09)       [0] [81] bestval =  6.45865e-09
##  1.000               0.004080                    6.45891e-09    (reg =            1            1   contrast =          6.45891e-09)       [0] [81] bestval =  6.45865e-09
->     3             6.45865e-09 <-          2.11015e-08 
~~~

Runtime 00:34:00 (size 1024)


---


## STEP 009 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput


This step takes approximately 4mn.


---


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


---


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
