/* zzgfspu.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static integer c__2 = 2;

/* $Procedure ZZGFSPU ( GF, angular separation utility routines ) */
/* Subroutine */ int zzgfspu_0_(int n__, char *of, char *from, char *shape, 
	char *frame, doublereal *refval, doublereal *et, char *abcorr, 
	logical *decres, logical *lssthn, doublereal *sep, ftnlen of_len, 
	ftnlen from_len, ftnlen shape_len, ftnlen frame_len, ftnlen 
	abcorr_len)
{
    /* Initialized data */

    static char svshap[32*2] = "POINT                           " "SPHERE   "
	    "                       ";
    static char ref[5] = "J2000";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int zzgftreb_(integer *, doublereal *);
    doublereal axes1[3], axes2[3];
    extern /* Subroutine */ int zzvalcor_(char *, logical *, ftnlen), chkin_(
	    char *, ftnlen), ucase_(char *, char *, ftnlen, ftnlen), errch_(
	    char *, char *, ftnlen, ftnlen);
    logical found;
    static doublereal svang;
    extern doublereal dvsep_(doublereal *, doublereal *);
    static char svref[32];
    static integer svobs;
    extern /* Subroutine */ int spkez_(integer *, doublereal *, char *, char *
	    , integer *, doublereal *, doublereal *, ftnlen, ftnlen), ljust_(
	    char *, char *, ftnlen, ftnlen), bods2c_(char *, integer *, 
	    logical *, ftnlen);
    static integer svbod1, svbod2;
    static doublereal svrad1, svrad2;
    static char svref1[32], svref2[32];
    extern logical failed_(void);
    static integer svshp1, svshp2;
    doublereal lt, dtheta;
    extern integer isrchc_(char *, integer *, char *, ftnlen, ftnlen);
    logical attblk[15];
    static char svabcr[32];
    extern doublereal zzdhfa_(doublereal *, doublereal *);
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen), cmprss_(char *, integer *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);
    doublereal seprtn;
    extern logical return_(void);
    doublereal pv1[6], pv2[6];
    extern /* Subroutine */ int zzgfspq_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, char *, char *, doublereal 
	    *, ftnlen, ftnlen);

/* $ Abstract */

/*     SPICE Private routine intended solely for the support of SPICE */
/*     routines. Users should not call this routine directly due */
/*     to the volatile nature of this routine. */

/*     This is the umbrella routine for the entry points needed by */
/*     GFEVNT in order to find angular separation events. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */
/* $ Abstract */

/*     Include file zzabcorr.inc */

/*     SPICE private file intended solely for the support of SPICE */
/*     routines.  Users should not include this file directly due */
/*     to the volatile nature of this file */

/*     The parameters below define the structure of an aberration */
/*     correction attribute block. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Parameters */

/*     An aberration correction attribute block is an array of logical */
/*     flags indicating the attributes of the aberration correction */
/*     specified by an aberration correction string.  The attributes */
/*     are: */

/*        - Is the correction "geometric"? */

/*        - Is light time correction indicated? */

/*        - Is stellar aberration correction indicated? */

/*        - Is the light time correction of the "converged */
/*          Newtonian" variety? */

/*        - Is the correction for the transmission case? */

/*        - Is the correction relativistic? */

/*    The parameters defining the structure of the block are as */
/*    follows: */

/*       NABCOR    Number of aberration correction choices. */

/*       ABATSZ    Number of elements in the aberration correction */
/*                 block. */

/*       GEOIDX    Index in block of geometric correction flag. */

/*       LTIDX     Index of light time flag. */

/*       STLIDX    Index of stellar aberration flag. */

/*       CNVIDX    Index of converged Newtonian flag. */

/*       XMTIDX    Index of transmission flag. */

/*       RELIDX    Index of relativistic flag. */

/*    The following parameter is not required to define the block */
/*    structure, but it is convenient to include it here: */

/*       CORLEN    The maximum string length required by any aberration */
/*                 correction string */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 18-DEC-2004 (NJB) */

/* -& */
/*     Number of aberration correction choices: */


/*     Aberration correction attribute block size */
/*     (number of aberration correction attributes): */


/*     Indices of attributes within an aberration correction */
/*     attribute block: */


/*     Maximum length of an aberration correction string: */


/*     End of include file zzabcorr.inc */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     OF         I   Names of the two targets */
/*     FROM       I   Name of the observing body */
/*     SHAPE      I   Names of the shape descriptions for OF */
/*     REFVAL     I   Anglular reference value for comparison */
/*     ET         I   An epoch in ephemeris seconds past J2000 TDB */
/*     ABCORR     I   Aberration correction flag */
/*     DECRES     O   .TRUE. if angular separation is decreasing .FALSE. */
/*                    otherwise */
/*     LSSTHN     O   .TRUE. is angular separation is less than REFVAL, */
/*                    .FALSE. otherwise */
/*     SEP        O   Angular separation at time ET */

/* $ Detailed_Input */

/*     OF       the string array naming the bodies whose angular */
/*              separation is of interest. */

/*     FROM     the string naming the observer. */

/*     SHAPE    the string array naming the geometric model used to */
/*              represent the shapes of OF. The relation between SHAPE */
/*              and OF is 1-to-1. */

/*              Models supported by this routine: */

/*                 'ELLIPSOID'     Use a triaxial ellipsoid model, */
/*                                 with radius values provided via the */
/*                                 kernel pool. A kernel variable */
/*                                 having a name of the form */

/*                                    'BODYnnn_RADII' */

/*                                 where nnn represents the NAIF */
/*                                 integer code associated with the */
/*                                 body, must be present in the kernel */
/*                                 pool. This variable must be */
/*                                 associated with three numeric */
/*                                 values giving the lengths of the */
/*                                 ellipsoid's X, Y, and Z semi-axes. */

/*                                 *This option not yet implemented.* */

/*                 'SPHERE'        Treat the body as a sphere with */
/*                                 radius equal to the maximum value of */
/*                                 BODYnnn_RADII */

/*                 'POINT'         Treat the body as a single point; */
/*                                 radius has value zero. */

/*              The SHAPE string lacks sensitivity to case and leading */
/*              or trailing blank. */

/*     FRAME    the string array naming the body-fixed reference frames */
/*              corresponding to OF. The relation between FRAME */
/*              and OF is 1-to-1. */

/*     REFVAL   the double precision value of the angle (in radians) */
/*              against which to compare the angular separation of the */
/*              two bodies. */

/*     ET       is the time in second past J2000 at which one wants */
/*              to determine an event condition. */

/*     ABCORR   the string description of the aberration corrections */
/*              to apply to the state evaluations to account for */
/*              one-way light time and stellar aberration. */

/*              This routine accepts the same aberration corrections */
/*              as does the SPICE routine SPKEZR. See the header of */
/*              SPKEZR for a detailed description of the aberration */
/*              correction options. For convenience, the options are */
/*              listed below: */

/*                 'NONE'     Apply no correction. */

/*                 'LT'       "Reception" case:  correct for */
/*                            one-way light time using a Newtonian */
/*                            formulation. */

/*                 'LT+S'     "Reception" case:  correct for */
/*                            one-way light time and stellar */
/*                            aberration using a Newtonian */
/*                            formulation. */

/*                 'CN'       "Reception" case:  converged */
/*                            Newtonian light time correction. */

/*                 'CN+S'     "Reception" case:  converged */
/*                            Newtonian light time and stellar */
/*                            aberration corrections. */

/*                 'XLT'      "Transmission" case:  correct for */
/*                            one-way light time using a Newtonian */
/*                            formulation. */

/*                 'XLT+S'    "Transmission" case:  correct for */
/*                            one-way light time and stellar */
/*                            aberration using a Newtonian */
/*                            formulation. */

/*                 'XCN'      "Transmission" case:  converged */
/*                            Newtonian light time correction. */

/*                 'XCN+S'    "Transmission" case:  converged */
/*                            Newtonian light time and stellar */
/*                            aberration corrections. */

/*                 The ABCORR string lacks sensitivity to case, leading */
/*                 and trailing blanks. */

/*     DECRES   is .TRUE. if the angular separation between the */
/*              objects is decreasing.  Otherwise it is .FALSE. */

/*     LSSTHN   is .TRUE. if the angular separation between the two */
/*              bodies is less than the reference angle at time ET */
/*              and .FALSE. otherwise. */

/*     SEP      is the angular separation between SVBOD1 and SVBOD2 as */
/*              seen from SVOBS at time ET. */

/*     For more information, see individual entry points. */

/* $ Detailed_Output */

/*     See individual entry points. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine serves as the umbrella routine for 4 entry points */
/*     needed by GFEVNT in solving for angular separation conditions. */

/*     The five entry points are */

/*        ZZGFSPIN  --- an initialization routine that must be called */
/*                     prior to attempting to solve for any angular */
/*                     separation event. */

/*        ZZGFSPUR --- updates reference value REFVAL. */

/*        ZZGFSPDC --- determines whether or not angular separation is */
/*                     decreasing at some time. */

/*        ZZGFSPLT --- determines whether or not angular separation is */
/*                     less than REFVAL */

/*        ZZGFGSEP --- returns the angular separation of the two */
/*                     objects of interest as a function of ET. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     ZZGFSPIN must be called prior to use of any of the */
/*     other entry points. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0 19-FEB-2009 (NJB) (EDW) */

/* -& */
/* $ Index_Entries */

/*     umbrella routine for finding angular separation events */

/* -& */

/*     SPICELIB functions */


/*     Local Variables */


/*     Saved Variables */


/*     Below we initialize the list of shape names. */


/*     Define integer ID parameters for the shape names in */
/*     SVSHAP. */

    /* Parameter adjustments */
    if (of) {
	}
    if (shape) {
	}
    if (frame) {
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_zzgfspin;
	case 2: goto L_zzgfspur;
	case 3: goto L_zzgfspdc;
	case 4: goto L_zzgfgsep;
	case 5: goto L_zzgfsplt;
	}


/*     Never directly call this routine. */

    chkin_("ZZGFSPU", (ftnlen)7);
    sigerr_("SPICE(BOGUSENTRY)", (ftnlen)17);
    chkout_("ZZGFSPU", (ftnlen)7);
    return 0;
/* $Procedure ZZGFSPIN (GF, angular separation initialization routine) */

L_zzgfspin:
/* $ Abstract */

/*     This routine initializes variables that describe an angular */
/*     separation event of interest for solution by ZZGFSOLV. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */

/*      CHARACTER*(*)         OF   ( 2 ) */
/*      INTEGER               FROM */
/*      CHARACTER*(*)         SHAPE( 2 ) */
/*      CHARACTER*(*)         FRAME( 2 ) */
/*      DOUBLE PRECISION      REFVAL */
/*      CHARACTER*(*)         ABCORR */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     OF         I   Body id's of the angular separation objects */
/*     FROM       I   Observer name */
/*     SHAPE      I   Array of shape IDs corresponding to OF */
/*     FRAME      I   Array of frame names corresponding to OF */
/*     REFVAL     I   Value angles will be compared to. */
/*     ABCORR     I   Aberration correction flag. */

/* $ Detailed_Input */

/*     OF       the string array naming the bodies whose angular */
/*              separation is of interest. */

/*     FROM     the string naming the observer. */

/*     SHAPE    the string array naming the geometric model used to */
/*              represent the shapes of OF. The relation between SHAPE */
/*              and OF is 1-to-1. */

/*              Models supported by this routine: */

/*                 'ELLIPSOID'     Use a triaxial ellipsoid model, */
/*                                 with radius values provided via the */
/*                                 kernel pool. A kernel variable */
/*                                 having a name of the form */

/*                                    'BODYnnn_RADII' */

/*                                 where nnn represents the NAIF */
/*                                 integer code associated with the */
/*                                 body, must be present in the kernel */
/*                                 pool. This variable must be */
/*                                 associated with three numeric */
/*                                 values giving the lengths of the */
/*                                 ellipsoid's X, Y, and Z semi-axes. */

/*                                 *This option not yet implemented.* */

/*                 'SPHERE'        Treat the body as a sphere with */
/*                                 radius equal to the maximum value of */
/*                                 BODYnnn_RADII */

/*                 'POINT'         Treat the body as a single point; */
/*                                 radius has value zero. */

/*              The SHAPE string lacks sensitivity to case and leading */
/*              or trailing blank. */

/*     FRAME    the string array naming the body-fixed reference frames */
/*              corresponding to OF. The relation between FRAME */
/*              and OF is 1-to-1. */

/*     REFVAL   the double precision value of the angle (in radians) */
/*              against which to compare the angular separation of the */
/*              two bodies. */

/*     ABCORR   the string description of the aberration corrections */
/*              to apply to the state evaluations to account for */
/*              one-way light time and stellar aberration. */

/*              This routine accepts the same aberration corrections */
/*              as does the SPICE routine SPKEZR. See the header of */
/*              SPKEZR for a detailed description of the aberration */
/*              correction options. For convenience, the options are */
/*              listed below: */

/*                 'NONE'     Apply no correction. */

/*                 'LT'       "Reception" case:  correct for */
/*                            one-way light time using a Newtonian */
/*                            formulation. */

/*                 'LT+S'     "Reception" case:  correct for */
/*                            one-way light time and stellar */
/*                            aberration using a Newtonian */
/*                            formulation. */

/*                 'CN'       "Reception" case:  converged */
/*                            Newtonian light time correction. */

/*                 'CN+S'     "Reception" case:  converged */
/*                            Newtonian light time and stellar */
/*                            aberration corrections. */

/*                 'XLT'      "Transmission" case:  correct for */
/*                            one-way light time using a Newtonian */
/*                            formulation. */

/*                 'XLT+S'    "Transmission" case:  correct for */
/*                            one-way light time and stellar */
/*                            aberration using a Newtonian */
/*                            formulation. */

/*                 'XCN'      "Transmission" case:  converged */
/*                            Newtonian light time correction. */

/*                 'XCN+S'    "Transmission" case:  converged */
/*                            Newtonian light time and stellar */
/*                            aberration corrections. */

/*                 The ABCORR string lacks sensitivity to case, leading */
/*                 and trailing blanks. */

/* $ Detailed_Output */

/*     None */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     None. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0 14-APR-2008 (NJB) (EDW) */

/* -& */
/* $ Index_Entries */

/*     angular separation initialization routine */

/* -& */
    if (return_()) {
	return 0;
    } else {
	chkin_("ZZGFSPIN", (ftnlen)8);
    }
    bods2c_(of, &svbod1, &found, of_len);
    if (! found) {
	setmsg_("The object name for target 1, '#', is not a recognized name"
		" for an ephemeris object. The cause of this problem may be t"
		"hat you need an updated version of the SPICE Toolkit.", (
		ftnlen)172);
	errch_("#", of, (ftnlen)1, of_len);
	sigerr_("SPICE(IDCODENOTFOUND)", (ftnlen)21);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }
    bods2c_(of + of_len, &svbod2, &found, of_len);
    if (! found) {
	setmsg_("The object name for target 2, '#', is not a recognized name"
		" for an ephemeris object. The cause of this problem may be t"
		"hat you need an updated version of the SPICE Toolkit.", (
		ftnlen)172);
	errch_("#", of + of_len, (ftnlen)1, of_len);
	sigerr_("SPICE(IDCODENOTFOUND)", (ftnlen)21);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }
    bods2c_(from, &svobs, &found, from_len);
    if (! found) {
	setmsg_("The object name for the observer, '#', is not a recognized "
		"name for an ephemeris object. The cause of this problem may "
		"be that you need an updated version of the SPICE Toolkit.", (
		ftnlen)176);
	errch_("#", from, (ftnlen)1, from_len);
	sigerr_("SPICE(IDCODENOTFOUND)", (ftnlen)21);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }

/*     Confirm the three bodies have unique IDs. */

    if (svobs == svbod1 || svobs == svbod2 || svbod1 == svbod2) {
	setmsg_("All three objects associated with an ANGULAR SEPARATION sea"
		"rch must be distinct. The objects whose angular separation i"
		"s of interest were # and #. The observer was #.", (ftnlen)166)
		;
	errint_("#", &svbod1, (ftnlen)1);
	errint_("#", &svbod2, (ftnlen)1);
	errint_("#", &svobs, (ftnlen)1);
	sigerr_("SPICE(BODIESNOTDISTINCT)", (ftnlen)24);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }

/*     Squeeze all blanks out of the aberration correction */
/*     string; ensure the string is in upper case. */

    cmprss_(" ", &c__0, abcorr, svabcr, (ftnlen)1, abcorr_len, (ftnlen)32);
    ucase_(svabcr, svabcr, (ftnlen)32, (ftnlen)32);

/*     Check the aberration correction. If SPKEZR can't handle it, */
/*     neither can we. */

    zzvalcor_(svabcr, attblk, (ftnlen)32);
    if (failed_()) {
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }
    s_copy(svref, ref, (ftnlen)32, (ftnlen)5);
    svang = *refval;
    s_copy(svref1, frame, (ftnlen)32, frame_len);
    s_copy(svref2, frame + frame_len, (ftnlen)32, frame_len);
    ljust_(shape, shape, shape_len, shape_len);
    ucase_(shape, shape, shape_len, shape_len);

/*     If we pass the error check, then SHAPE(1) */
/*     exists in SVSHAP. */

    svshp1 = isrchc_(shape, &c__2, svshap, shape_len, (ftnlen)32);
    if (svshp1 == 0) {
	setmsg_("The body shape, # is not recognized.  Supported quantities "
		"are: POINT, SPHERE.", (ftnlen)78);
	errch_("#", shape, (ftnlen)1, shape_len);
	sigerr_("SPICE(NOTRECOGNIZED)", (ftnlen)20);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    } else if (svshp1 == 1) {
	svrad1 = 0.;
    } else if (svshp1 == 2) {
	zzgftreb_(&svbod1, axes1);
	if (failed_()) {
	    chkout_("ZZGFSPIN", (ftnlen)8);
	    return 0;
	}
/* Computing MAX */
	d__1 = max(axes1[0],axes1[1]);
	svrad1 = max(d__1,axes1[2]);
    } else {

/*        This code executes only if someone adds a new shape */
/*        name to SVSHAP then fails to update the SVSHP1 condition */
/*        block to respond to the name. Fortran needs SWITCH...CASE. */

	setmsg_("Encountered uncoded shape ID for #. This indicates a bog. P"
		"lease contact NAIF.", (ftnlen)78);
	errch_("#", shape, (ftnlen)1, shape_len);
	sigerr_("SPICE(BUG)", (ftnlen)10);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }
    ljust_(shape + shape_len, shape + shape_len, shape_len, shape_len);
    ucase_(shape + shape_len, shape + shape_len, shape_len, shape_len);

/*     If we pass the error check, then SHAPE(2) */
/*     exists in SVSHAP. */

    svshp2 = isrchc_(shape + shape_len, &c__2, svshap, shape_len, (ftnlen)32);
    if (svshp2 == 0) {
	setmsg_("The body shape, # is not recognized.  Supported quantities "
		"are: POINT, SPHERE.", (ftnlen)78);
	errch_("#", shape + shape_len, (ftnlen)1, shape_len);
	sigerr_("SPICE(NOTRECOGNIZED)", (ftnlen)20);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    } else if (svshp2 == 1) {
	svrad2 = 0.;
    } else if (svshp2 == 2) {
	zzgftreb_(&svbod2, axes2);
	if (failed_()) {
	    chkout_("ZZGFSPIN", (ftnlen)8);
	    return 0;
	}
/* Computing MAX */
	d__1 = max(axes2[0],axes2[1]);
	svrad2 = max(d__1,axes2[2]);
    } else {

/*        This code executes only if someone adds a new shape */
/*        name to SVSHAP then fails to update the SVSHP2 condition */
/*        block to respond to the name. Fortran needs SWITCH...CASE. */

	setmsg_("Encountered uncoded shape ID for #. This indicates a bug. P"
		"lease contact NAIF.", (ftnlen)78);
	errch_("#", shape + shape_len, (ftnlen)1, shape_len);
	sigerr_("SPICE(BUG)", (ftnlen)10);
	chkout_("ZZGFSPIN", (ftnlen)8);
	return 0;
    }
    chkout_("ZZGFSPIN", (ftnlen)8);
    return 0;
/* $Procedure ZZGFSPUR (GF, update angular separation reference value ) */

L_zzgfspur:
/* $ Abstract */

/*     This is the entry point used for updating the internal reference */
/*     value. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */

/*      DOUBLE PRECISION      REFVAL */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     REFVAL     I   Anglular reference value for comparison */

/* $ Detailed_Input */

/*     REFVAL     the double precision value of the angle (in radians) */
/*                against which to compare the angular separation of the */
/*                two bodies. */

/* $ Detailed_Output */

/*     None */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     None. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 17-FEB-2009 (EDW) */

/* -& */
/* $ Index_Entries */

/*     angular separation update reference value routine */

/* -& */
    svang = *refval;
    return 0;
/* $Procedure ZZGFSPDC (GF, angular separation decreasing) */

L_zzgfspdc:
/* $ Abstract */

/*     Computes whether or not the angular separation between SVBOD1 and */
/*     SVBOD2 is decreasing at time ET. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */

/*     DOUBLE PRECISION      ET */
/*     LOGICAL               DECRES */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     ET         I   Ephemeris seconds past J2000 TDB. */
/*     DECRES     O   .TRUE if angular separation is decreasing .FALSE. */
/*                    otherwise. */

/* $ Detailed_Input */

/*     ET         time in seconds past J2000 at which one wishes to */
/*                determine whether or not the angular separation of the */
/*                two bodies is decreasing. */

/* $ Detailed_Output */

/*     DECRES     is .TRUE. if the angular separation between the objects */
/*                is decreasing.  Otherwise it is .FALSE. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     If the observer is inside one of the objects, the object will */
/*     be regarded as having a 90 degree apparent radius. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine determines whether or not the angular separation */
/*     between two objects as seen from a third is decreasing. The value */
/*     of DECRES is .TRUE. if it is, otherwise it is returned as */
/*     .FALSE. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */
/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0 29-APR-2008 (NJB) */

/* -& */
/* $ Index_Entries */

/*     angular separation is decreasing */

/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("ZZGFSPDC", (ftnlen)8);
    }
    spkez_(&svbod1, et, svref, svabcr, &svobs, pv1, &lt, (ftnlen)32, (ftnlen)
	    32);
    if (failed_()) {
	chkout_("ZZGFSPDC", (ftnlen)8);
	return 0;
    }
    spkez_(&svbod2, et, svref, svabcr, &svobs, pv2, &lt, (ftnlen)32, (ftnlen)
	    32);
    if (failed_()) {
	chkout_("ZZGFSPDC", (ftnlen)8);
	return 0;
    }

/*     The angular separation between the bodies has the value */

/*        theta = sep - alpha1 - alpha2 */

/*     With alpha1 the half angle of SVBOD1, alpha2 the half */
/*     angle of SVBOD2, half angle defined as (for spheres): */

/*        sin(alpha) = body_radius */
/*                     ----------- */
/*                     range_to_body */

/*     The corresponding time derivative of theta: */

/*        d(theta) = d(sep) - d(alpha1) - d(alpha2) */
/*        --------   ------   ---------   --------- */
/*        dt         dt       dt          dt */

/*     Note, alpha1, alpha2 and their derivatives have value zero */
/*     for point objects. */

    dtheta = dvsep_(pv1, pv2);

/*     Check for a failure caused by a numerical event. */

    if (failed_()) {
	*decres = TRUE_;
	chkout_("ZZGFSPDC", (ftnlen)8);
	return 0;
    }
    dtheta = dtheta - zzdhfa_(pv1, &svrad1) - zzdhfa_(pv2, &svrad2);
    if (dtheta < 0.) {
	*decres = TRUE_;
    } else {
	*decres = FALSE_;
    }
    chkout_("ZZGFSPDC", (ftnlen)8);
    return 0;
/* $Procedure ZZGFGSEP (GF, calculate angular separation between bodies ) */

L_zzgfgsep:
/* $ Abstract */

/*     Determine the angular separation between the limbs of the two */
/*     bodies. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */

/*      DOUBLE PRECISION      ET */
/*      DOUBLE PRECISION      SEP */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     ET         I   Ephemeris seconds past J2000 TDB. */
/*     SEP        O   Separation at time ET. */

/* $ Detailed_Input */

/*     ET         time in ephemeris seconds past J2000 when the */
/*                angular separation between the two bodies is */
/*                to be computed. */

/* $ Detailed_Output */

/*     SEP        is the angular separation between SVBOD1 and SVBOD2 as */
/*                seen from SVOBS at time ET. */


/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine determins the apparent angular separation between the */
/*     limbs of bodies SVBOD1 and SVBOD2 as seen from SVOBS at time ET. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0 26-AUG-2003 (LSE) */

/* -& */
/* $ Index_Entries */

/*     angular separation between two bodies */

/* -& */
    zzgfspq_(et, &svbod1, &svbod2, &svrad1, &svrad2, &svobs, svabcr, svref, 
	    sep, (ftnlen)32, (ftnlen)32);
    return 0;
/* $Procedure ZZGFSPLT  (GF, angular separation less than reference ) */

L_zzgfsplt:
/* $ Abstract */

/*     Determine whether or not the angular separation between the two */
/*     bodies is less than the reference value. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     ANGLE */
/*     GEOMETRY */
/*     ROOT */

/* $ Declarations */

/*      DOUBLE PRECISION      ET */
/*      LOGICAL               LSSTHN */

/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     ET         I   Ephemeris seconds past J2000 TDB. */
/*     LSSTHN     O   .TRUE. if separation is less than REFVAL, */
/*                    .FALSE. otherwise. */

/* $ Detailed_Input */

/*     ET         is the time in second past J2000 at which one wants */
/*                to determine if the angular separation between the */
/*                two bodies is less than the reference angle. */

/* $ Detailed_Output */

/*     LSSTHN     is .TRUE. if the angle between the two bodies is less */
/*                than the reference angle at time ET and .FALSE. */
/*                otherwise. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     None. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This routine determines whether or not the angle between */
/*     the two objects as seen from SVOBS is less than the reference */
/*     angle at time ET. */

/* $ Examples */

/*     None. */

/* $ Restrictions */

/*     1) Due to the current logic implemented in ZZGFSPU, a direct */
/*        search for the zero angular separation of two point targets */
/*        will always fails, i.e., */

/*           OP     = '=' */
/*           REFVAL = 0.D0. */

/*        Use OP values of 'ABSMIN' or 'LOCMIN' to detect such an event. */


/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     W.L. Taber     (JPL) */
/*     I.M. Underwood (JPL) */
/*     L.S. Elson     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0 19-FEB-2009 (EDW) */

/* -& */
/* $ Index_Entries */

/*     angular separation less than an angle */

/* -& */
    zzgfspq_(et, &svbod1, &svbod2, &svrad1, &svrad2, &svobs, svabcr, svref, &
	    seprtn, (ftnlen)32, (ftnlen)32);
    if (seprtn < svang) {
	*lssthn = TRUE_;
    } else {
	*lssthn = FALSE_;
    }
    return 0;
} /* zzgfspu_ */

/* Subroutine */ int zzgfspu_(char *of, char *from, char *shape, char *frame, 
	doublereal *refval, doublereal *et, char *abcorr, logical *decres, 
	logical *lssthn, doublereal *sep, ftnlen of_len, ftnlen from_len, 
	ftnlen shape_len, ftnlen frame_len, ftnlen abcorr_len)
{
    return zzgfspu_0_(0, of, from, shape, frame, refval, et, abcorr, decres, 
	    lssthn, sep, of_len, from_len, shape_len, frame_len, abcorr_len);
    }

/* Subroutine */ int zzgfspin_(char *of, char *from, char *shape, char *frame,
	 doublereal *refval, char *abcorr, ftnlen of_len, ftnlen from_len, 
	ftnlen shape_len, ftnlen frame_len, ftnlen abcorr_len)
{
    return zzgfspu_0_(1, of, from, shape, frame, refval, (doublereal *)0, 
	    abcorr, (logical *)0, (logical *)0, (doublereal *)0, of_len, 
	    from_len, shape_len, frame_len, abcorr_len);
    }

/* Subroutine */ int zzgfspur_(doublereal *refval)
{
    return zzgfspu_0_(2, (char *)0, (char *)0, (char *)0, (char *)0, refval, (
	    doublereal *)0, (char *)0, (logical *)0, (logical *)0, (
	    doublereal *)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (
	    ftnint)0);
    }

/* Subroutine */ int zzgfspdc_(doublereal *et, logical *decres)
{
    return zzgfspu_0_(3, (char *)0, (char *)0, (char *)0, (char *)0, (
	    doublereal *)0, et, (char *)0, decres, (logical *)0, (doublereal *
	    )0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0);
    }

/* Subroutine */ int zzgfgsep_(doublereal *et, doublereal *sep)
{
    return zzgfspu_0_(4, (char *)0, (char *)0, (char *)0, (char *)0, (
	    doublereal *)0, et, (char *)0, (logical *)0, (logical *)0, sep, (
	    ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0);
    }

/* Subroutine */ int zzgfsplt_(doublereal *et, logical *lssthn)
{
    return zzgfspu_0_(5, (char *)0, (char *)0, (char *)0, (char *)0, (
	    doublereal *)0, et, (char *)0, (logical *)0, lssthn, (doublereal *
	    )0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0, (ftnint)0);
    }

