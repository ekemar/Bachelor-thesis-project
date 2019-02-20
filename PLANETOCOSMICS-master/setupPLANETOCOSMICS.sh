#! /bin/bash 
#Setup script for the PLANETOCOSMICS code
#Written by L. Desorgher
export IGRF_TABLE=${PLANETOCOSMICS}/data/earth/igrfdata/igrf_dgrf_table.txt
export PATH=${PLANETOCOSMICS}/bin:$PATH:
export MARS_CRUSTAL=${PLANETOCOSMICS}/data/mars/fieldtables
export SPICE_LEAPSECOND_DATA=${PLANETOCOSMICS}/lib/cspice/data/leap_seconds.tls
export SPICE_FRAME_DATA=${PLANETOCOSMICS}/lib/cspice/data/pck00008.tpc
export SPICE_EPHEMERID_DATA=${PLANETOCOSMICS}/lib/cspice/data/de410.bsp
export MARS_ATMO_DATA=${PLANETOCOSMICS}/data/mars/marsgram_atmotables
