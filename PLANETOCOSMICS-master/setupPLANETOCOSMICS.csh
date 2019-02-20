#! /bin/tcsh 
#Setup script for the version of PLANETOCOSMICS that uses ROOT
#Written by L. Desorgher
setenv IGRF_TABLE ${PLANETOCOSMICS}/data/earth/igrfdata/igrf_dgrf_table.txt
setenv PATH  ${PLANETOCOSMICS}/bin:${PATH} 
setenv MARS_CRUSTAL ${PLANETOCOSMICS}/data/mars/fieldtables 
setenv SPICE_LEAPSECOND_DATA ${PLANETOCOSMICS}/lib/cspice/data/leap_seconds.tls
setenv SPICE_FRAME_DATA ${PLANETOCOSMICS}/lib/cspice/data/pck00008.tpc
setenv SPICE_EPHEMERID_DATA ${PLANETOCOSMICS}/lib/cspice/data/de410.bsp
setenv MARS_ATMO_DATA ${PLANETOCOSMICS}/data/mars/marsgram_atmotables 
 
