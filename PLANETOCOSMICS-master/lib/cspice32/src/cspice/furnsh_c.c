/*

-Procedure furnsh_c ( Furnish a program with SPICE kernels )

-Abstract
 
   Load one or more SPICE kernels into a program.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading
 
    None. 
 
-Keywords
 
    UTILITY 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"


   void furnsh_c ( ConstSpiceChar  * file ) 

/*

-Brief_I/O
 
   VARIABLE  I/O  DESCRIPTION 
   --------  ---  -------------------------------------------------- 
   file       I   Name of SPICE kernel file (text or binary). 
 
-Detailed_Input
 
   file       is the name of a SPICE kernel file. The file may be
              either binary or text. If the file is a binary SPICE
              kernel it will be loaded into the appropriate SPICE
              subsystem.  If `file' is a SPICE text kernel it will be
              loaded into the kernel pool.  If `file' is a SPICE
              meta-kernel containing initialization instructions
              (through use of the correct kernel pool variables), the
              files specified in those variables will be loaded into
              the appropriate SPICE subsystem.

              The SPICE text kernel format supports association of
              names and data values using a "keyword = value" format.
              The keyword-value pairs thus defined are called "kernel
              variables." 

              While any information can be placed in a text kernel 
              file, the following string valued kernel variables are 
              recognized by SPICE as meta-kernel keywords:
 
                 KERNELS_TO_LOAD 
                 PATH_SYMBOLS 
                 PATH_VALUES 
 
              Each kernel variable is discussed below. 
 
                 KERNELS_TO_LOAD   is a list of SPICE kernels to be 
                                   loaded into a program. If file 
                                   names do not fit within the kernel 
                                   pool 80 character limit, they may be 
                                   continued to subsequent array 
                                   elements by placing the continuation 
                                   character ('+') at the end of an 
                                   element and then placing the 
                                   remainder of the file name in the 
                                   next array element. (See the 
                                   examples below for an illustration 
                                   of this technique or consult the 
                                   routine stpool_c for further 
                                   details.) 
    
                                   You may use one or more PATH_SYMBOL
                                   assignments (see below) to specify
                                   strings to be substituted for some
                                   part of a file name.
    
                 PATH_SYMBOLS      is a list of strings (without 
                                   embedded blanks), which if 
                                   encountered following the '$' 
                                   character will be replaced with the 
                                   corresponding PATH_VALUES string. 
                                   Note that PATH_SYMBOLS are
                                   interpreted only in values
                                   associated with the KERNELS_TO_LOAD
                                   variable. There must be a one-to-one
                                   correspondence between the values
                                   supplied for PATH_SYMBOLS and
                                   PATH_VALUES. For the purpose of
                                   determining this correspondence, any
                                   path value that is continued over
                                   multiple array elements counts as a
                                   single value.
    
                 PATH_VALUES       is a list of expansions to use when 
                                   PATH_SYMBOLS are encountered. If
                                   path values do not fit within the
                                   kernel pool 80 character limit, they
                                   may be continued in the same way as
                                   file names (see the KERNELS_TO_LOAD
                                   description above).

 
              These kernel pool variables persist within the kernel 
              pool only until all kernels associated with the 
              variable KERNELS_TO_LOAD have been loaded.  Once all 
              specified kernels have been loaded, the variables 
              KERNELS_TO_LOAD, PATH_SYMBOLS and PATH_VALUES are 
              removed from the kernel pool. 
 
-Detailed_Output
 
   None. The routine loads various SPICE kernels for use by your 
   application. 
 
-Parameters
 
   FILSIZ     is the maximum length of a file name that can
              be loaded by this routine. FILSIZ is currently
              set to 255 characters.
  
-Exceptions
 

   1) If a problem is encountered while trying to load `file',
      it will be diagnosed by a routine in the call tree of this
      routine.

   2) If the input `file' is a meta-kernel and some file in the
      KERNELS_TO_LOAD assignment cannot be found, or if an error
      occurs while trying to load a file specified by this
      assignment, the error will be diagnosed by a routine in the
      call tree of this routine, and this routine will return. Any
      files loaded prior to encountering the missing file will
      remain loaded.

   3) If a PATH_SYMBOLS assignment is specified without a
      corresponding PATH_VALUES assignment, the error
      SPICE(NOPATHVALUE) will be signaled.

   4) If a meta-text kernel is supplied to furnsh_c that contains
      instructions specifying that another meta-text kernel be
      loaded, the error SPICE(RECURSIVELOADING) will be signaled.

   5) If the input file name has non-blank length exceeding FILSIZ
      characters, the error SPICE(FILENAMETOOLONG) is signaled.

   6) If the input file is a meta-kernel and some file in the
      KERNELS_TO_LOAD assignment has name length exceeding FILSIZ
      characters, the error SPICE(FILENAMETOOLONG) is signaled.

   7) If the input file is a meta-kernel and some value in the
      PATH_VALUES assignment has length exceeding FILSIZ
      characters, the error SPICE(PATHTOOLONG) is signaled.

   8) If the input file is a meta-kernel and some file in the
      KERNELS_TO_LOAD assignment has, after symbol substitution,
      combined name and path length exceeding FILSIZ characters, the
      error SPICE(FILENAMETOOLONG) is signaled.

   9) If the input `file' argument pointer is null, the error
      SPICE(NULLPOINTER) will be signaled.
      
   10) If the input `file' argument is the empty string, the error
       SPICE(EMPTYSTRING) will be signaled.

   
-Files
 
   The input file is examined and loaded into the appropriate 
   SPICE subsystem.  If the file is a meta-kernel, any kernels
   specified by the KERNELS_TO_LOAD keyword (and if present,
   the PATH_SYMBOLS and PATH_VALUES keywords) are loaded as well.
   
-Particulars
 
   This routine provides a uniform interface to the SPICE kernel 
   loading systems.  It allows you to easily assemble a list of 
   SPICE kernels required by your application and to modify that set 
   without modifying the source code of programs that make use of 
   these kernels. 

   Text kernels input to this routine need not have native line
   terminators for the platform. Lower level CSPICE routines can
   read and process non-native text files. This functionality does
   not exist in the Fortran SPICELIB.

   Only text kernel readers include the non-native read capability,
   (ldpool_c and furnsh_c), the generic text file line reader, rdtext_c
   requires native text files.
   
   Please refer to kernel.req for additional information.
    
-Examples
 
   Example 1
   --------- 

   Load the leapseconds kernel naif0007.tls and the planetary ephemeris
   SPK file de405s.bsp.

       furnsh_c ( "naif0007.tls" );
       furnsh_c ( "de405s.bsp"   );


   Example 2 
   --------- 
 
   This example illustrates how you could create a meta-kernel file for
   a program that requires several text and binary kernels.
 
   First create a list of the kernels you need in a text file as
   shown below.
 
      \begintext 
 
         Here are the SPICE kernels required for my application 
         program. 
 
         Note that kernels are loaded in the order listed. Thus we 
         need to list the highest priority kernel last. 
 
 
      \begindata 
      
      KERNELS_TO_LOAD = ( '/home/mydir/kernels/spk/lowest_priority.bsp', 
                          '/home/mydir/kernels/spk/next_priority.bsp', 
                          '/home/mydir/kernels/spk/highest_priority.bsp', 
                          '/home/mydir/kernels/text/leapsecond.ker', 
                          '/home/mydir/kernels+',
                          '/custom+', 
                          '/kernel_data/constants.ker',
                          '/home/mydir/kernels/text/sclk.tsc', 
                          '/home/mydir/kernels/ck/c-kernel.bc' ) 
 
 
   Note that the file name
 
      /home/mydir/kernels/custom/kernel_data/constants.ker 
 
   is continued across several lines in the right hand side of the
   assignment of the kernel variable KERNELS_TO_LOAD. 
 
   Once you've created your list of kernels, call furnsh_c near the
   beginning of your application program to load the meta-kernel
   automatically at program start up.
 
      furnsh_c ( "myfile.txt" );

   This will cause each of the kernels listed in your meta-kernel
   to be loaded.
 
 
   Example 3
   --------- 
 
   This example illustrates how you can simplify the previous 
   kernel list by using PATH_SYMBOLS. 
 
 
      \begintext 
 
         Here are the SPICE kernels required for my application 
         program. 
 
         We are going to let A substitute for the directory that
         contains SPK files; B substitute for the directory that
         contains C-kernels; and C substitute for the directory that
         contains text kernels. And we'll let D substitute for
         a "custom" directory that contains a special planetary
         constants kernel made just for our mission.
 
         Note that our PATH_VALUES and the corresponding
         PATH_SYMBOLS must be listed in the same order.
 
      \begindata 
      
      PATH_VALUES  = ( '/home/mydir/kernels/spk', 
                       '/home/mydir/kernels/ck', 
                       '/home/mydir/kernels/text',
                       '/home/mydir/kernels/custom/kernel_data' ) 
 
      PATH_SYMBOLS = ( 'A', 
                       'B', 
                       'C'
                       'D'  ) 
 
 
      KERNELS_TO_LOAD = (  '$A/lowest_priority.bsp', 
                           '$A/next_priority.bsp', 
                           '$A/highest_priority.bsp', 
                           '$C/leapsecond.ker', 
                           '$D/constants.ker',
                           '$C/sclk.tsc', 
                           '$B/c-kernel.bc'         )
 
 
   Example 4
   ---------

   This example illustrates continuation of path values. The
   meta-kernel shown here is a modified version of that from
   example 3.

      \begintext

         Here are the SPICE kernels required for my application
         program.

         We are going to let A substitute for the directory that
         contains SPK files; B substitute for the directory that
         contains C-kernels; and C substitute for the directory that
         contains text kernels. And we'll let D substitute for
         a "custom" directory that contains a special planetary
         constants kernel made just for our mission.

         Note that our PATH_VALUES and the corresponding
         PATH_SYMBOLS must be listed in the same order.

         The values for path symbols A and D are continued over
         multiple lines.

      \begindata

      PATH_VALUES  = ( '/very_long_top_level_path_name/mydir/+',
                       'kernels/spk',
                       '/home/mydir/kernels/ck',
                       '/home/mydir/kernels/text',
                       '/very_long_top_level_path_name+',
                       '/mydir/kernels/custom+',
                       '/kernel_data'                )

      PATH_SYMBOLS = ( 'A',
                       'B',
                       'C',
                       'D'  )

      KERNELS_TO_LOAD = (  '$A/lowest_priority.bsp',
                           '$A/next_priority.bsp',
                           '$A/highest_priority.bsp',
                           '$C/leapsecond.ker',
                           '$D/constants.ker',
                           '$C/sclk.tsc',
                           '$B/c-kernel.bc'         )


-Restrictions
 
   None. 
 
-Literature_References

   None. 
 
-Author_and_Institution
   
   C.H. Acton      (JPL)
   N.J. Bachman    (JPL)
   W.L. Taber      (JPL) 
   E.D. Wright     (JPL)
 
-Version

   -CSPICE Version 2.0.0, 02-APR-2009 (NJB)
 
        Continued path values are now supported. furnsh_c now rejects
        file names longer than FILSIZ characters.

   -CSPICE Version 1.0.4 17-OCT-2005 (EDW)

      Added text to Particulars section informing of the
      non-native kernel text file reading capability.

   -CSPICE Version 1.0.3, 29-JUL-2003 (NJB) (CHA)

      Numerous updates to improve clarity.  Some corrections
      were made.
 
   -CSPICE Version 1.0.2, 03-JUL-2002 (NJB) 
   
      Documentation fix:  corrected second code example.  The example
      previously used the kernel variable PATH_NAMES; that name has
      been replaced with the correct name PATH_VALUES.
 
   -CSPICE Version 1.0.1, 13-APR-2000 (NJB) 
   
      Replaced single quotes with double quotes in a code example.

   -CSPICE Version 1.0.0, 01-SEP-1999 (NJB) (WLT)

-Index_Entries
 
   Load SPICE data from a list of items 
 
-&
*/

{ /* Begin furnsh_c */



   /*
   Participate in error tracing.
   */
   chkin_c ( "furnsh_c" );


   /*
   Check the input filename to make sure the pointer is non-null 
   and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "furnsh_c", file );


   /*
   Call the f2c'd Fortran routine.
   */
   furnsh_ ( ( char   * ) file, 
             ( ftnlen   ) strlen(file) );


   chkout_c ( "furnsh_c" );

} /* End furnsh_c */
