#ifndef nrlmsise00h
#define nrlmsise00h

/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/* This file is part of the NRLMSISE-00  C source code package - release
 * 20041227
 *
 * The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
 * Doug Drob. They also wrote a NRLMSISE-00 distribution package in 
 * FORTRAN which is available at
 * http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
 *
 * Dominik Brodowski implemented and maintains this C version. You can
 * reach him at mail@brodo.de. See http://www.brodo.de/english/pub/nrlmsise/index.html for
 * updated releases of this package.
 *
 * C++ Version by Philip von Doetinchem, doetinchem@ssl.berkeley.edu
 */



/* ------------------------------------------------------------------- */
/* ------------------------------- INPUT ----------------------------- */
/* ------------------------------------------------------------------- */

class nrlmsiseflags 
{
public:
int switches[24];   
double sw[24];
double swc[24];
};

/*   
 *   Switches: to turn on and off particular variations use these switches.
 *   0 is off, 1 is on, and 2 is main effects off but cross terms on.
 *
 *   Standard values are 0 for switch 0 and 1 for switches 1 to 23. The 
 *   array "switches" needs to be set accordingly by the calling program. 
 *   The arrays sw and swc are set internally.
 *
 *   switches[i]:
 *    i - explanation
 *   -----------------
 *    0 - output in centimeters instead of meters
 *    1 - F10.7 effect on mean
 *    2 - time independent
 *    3 - symmetrical annual
 *    4 - symmetrical semiannual
 *    5 - asymmetrical annual
 *    6 - asymmetrical semiannual
 *    7 - diurnal
 *    8 - semidiurnal
 *    9 - daily ap [when this is set to -1 (!) the pointer
 *                  ap_a in struct nrlmsise_input must
 *                  point to a struct ap_array]
 *   10 - all UT/long effects
 *   11 - longitudinal
 *   12 - UT and mixed UT/long
 *   13 - mixed AP/UT/LONG
 *   14 - terdiurnal
 *   15 - departures from diffusive equilibrium
 *   16 - all TINF var
 *   17 - all TLB var
 *   18 - all TN1 var
 *   19 - all S var
 *   20 - all TN2 var
 *   21 - all NLB var
 *   22 - all TN3 var
 *   23 - turbo scale height var
 */

class aparray 
{
public:
double a[7];   
};
/* Array containing the following magnetic values:
 *   0 : daily AP
 *   1 : 3 hr AP index for current time
 *   2 : 3 hr AP index for 3 hrs before current time
 *   3 : 3 hr AP index for 6 hrs before current time
 *   4 : 3 hr AP index for 9 hrs before current time
 *   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
 *           prior to current time
 *   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
 *           prior to current time 
 */


class nrlmsiseinput 
{
public:
int year;      /* year, currently ignored */
int doy;       /* day of year */
double sec;    /* seconds in day (UT) */
double alt;    /* altitude in kilometes */
double g_lat;  /* geodetic latitude */
double g_long; /* geodetic longitude */
double lst;    /* local apparent solar time (hours), see note below */
double f107A;  /* 81 day average of F10.7 flux (centered on doy) */
double f107;   /* daily F10.7 flux for previous day */
double ap;     /* magnetic index(daily) */
aparray *ap_a; /* see above */
};
/*
 *   NOTES ON INPUT VARIABLES: 
 *      UT, Local Time, and Longitude are used independently in the
 *      model and are not of equal importance for every situation.  
 *      For the most physically realistic calculation these three
 *      variables should be consistent (lst=sec/3600 + g_long/15).
 *      The Equation of Time departures from the above formula
 *      for apparent local time can be included if available but
 *      are of minor importance.
 *
 *      f107 and f107A values used to generate the model correspond
 *      to the 10.7 cm radio flux at the actual distance of the Earth
 *      from the Sun rather than the radio flux at 1 AU. The following
 *      site provides both classes of values:
 *      ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
 *
 *      f107, f107A, and ap effects are neither large nor well
 *      established below 80 km and these parameters should be set to
 *      150., 150., and 4. respectively.
 */



/* ------------------------------------------------------------------- */
/* ------------------------------ OUTPUT ----------------------------- */
/* ------------------------------------------------------------------- */

class nrlmsiseoutput 
{
public:
double d[9];   /* densities */
double t[2];   /* temperatures */
};
/* 
 *   OUTPUT VARIABLES:
 *      d[0] - HE NUMBER DENSITY(CM-3)
 *      d[1] - O NUMBER DENSITY(CM-3)
 *      d[2] - N2 NUMBER DENSITY(CM-3)
 *      d[3] - O2 NUMBER DENSITY(CM-3)
 *      d[4] - AR NUMBER DENSITY(CM-3)                       
 *      d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
 *      d[6] - H NUMBER DENSITY(CM-3)
 *      d[7] - N NUMBER DENSITY(CM-3)
 *      d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
 *      t[0] - EXOSPHERIC TEMPERATURE
 *      t[1] - TEMPERATURE AT ALT
 * 
 *
 *      O, H, and N are set to zero below 72.5 km
 *
 *      t[0], Exospheric temperature, is set to global average for
 *      altitudes below 120 km. The 120 km gradient is left at global
 *      average value for altitudes below 72 km.
 *
 *      d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
 *      and GTD7D
 *
 *        SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
 *        species labeled by indices 0-4 and 6-7 in output variable d.
 *        This includes He, O, N2, O2, Ar, H, and N but does NOT include
 *        anomalous oxygen (species index 8).
 *
 *        SUBROUTINE GTD7D -- d[5] is the "effective total mass density
 *        for drag" and is the sum of the mass densities of all species
 *        in this model, INCLUDING anomalous oxygen.
 */



class nrlmsisecalc 
{
public:
nrlmsisecalc();
~nrlmsisecalc();
	
/* GTD7 */
/*   Neutral Atmosphere Empircial Model from the surface to lower
 *   exosphere.
 */
void gtd7 (nrlmsiseinput *input, nrlmsiseflags *flags, nrlmsiseoutput *output);

/* GTD7D */
/*   This subroutine provides Effective Total Mass Density for output
 *   d[5] which includes contributions from "anomalous oxygen" which can
 *   affect satellite drag above 500 km. See the section "output" for
 *   additional details.
 */
void gtd7d(nrlmsiseinput *input, nrlmsiseflags *flags,nrlmsiseoutput *output);
	   
/* GTS7 */
/*   Thermospheric portion of NRLMSISE-00
 */
void gts7 (nrlmsiseinput *input, nrlmsiseflags *flags, nrlmsiseoutput *output);
	   
/* GHP7 */
/*   To specify outputs at a pressure level (press) rather than at
 *   an altitude.
 */
void ghp7 (nrlmsiseinput *input, nrlmsiseflags *flags, nrlmsiseoutput *output, double press);

void tselec(nrlmsiseflags *flags);
void glatf(double lat, double *gv, double *reff);
double ccor(double alt, double r, double h1, double zh);
double ccor2(double alt, double r, double h1, double zh, double h2);
double scalh(double alt, double xm, double temp);
double dnet (double dd, double dm, double zhm, double xmm, double xm);
void splini (double *xa, double *ya, double *y2a, int n, double x, double *y);
void splint (double *xa, double *ya, double *y2a, int n, double x, double *y);
void spline (double *x, double *y, int n, double yp1, double ypn, double *y2);
double zeta(double zz, double zl);
double densm (double alt, double d0, double xm, double *tz, int mn3, double *zn3, double *tn3, double *tgn3, int mn2, double *zn2, double *tn2, double *tgn2);
double densu (double alt, double dlb, double tinf, double tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double *zn1, double *tn1, double *tgn1);
double g0(double a, double *p);
double sumex(double ex);
double sg0(double ex, double *p, double *ap);
double globe7(double *p, nrlmsiseinput *input, nrlmsiseflags *flags);
double glob7s(double *p, nrlmsiseinput *input, nrlmsiseflags *flags);

private:
/* PARMB */
double gsurf;
double re;

/* GTS3C */
double dd;

/* DMIX */
double dm04; 
double dm16; 
double dm28; 
double dm32; 
double dm40; 
double dm01;
double dm14;

/* MESO7 */
double meso_tn1[5];
double meso_tn2[4];
double meso_tn3[5];
double meso_tgn1[2];
double meso_tgn2[2];
double meso_tgn3[2];

/* POWER7 */
double pt[150];
double pd[9][150];
double ps[150];
double pdl[2][25];
double ptl[4][100];
double pma[10][100];
double sam[100];

/* LOWER7 */
double ptm[10];
double pdm[8][10];
double pavgm[10];

/* LPOLY */
double dfa;
double plg[4][9];
double ctloc; 
double stloc;
double c2tloc; 
double s2tloc;
double s3tloc; 
double c3tloc;
double apdf; 
double apt[4]; 
};

#endif












