#! /bin/sh

 ncap -O -s "TOTn[time, depth, lat, lon]=B1n+N3n+N4n+O4n+P1n+P2n+P3n+P4n+R1n+R6n+Z3n+Z4n+Z5n+Z6n; Nit[time, depth, lat, lon]=N3n+N4n; ORGn[time, depth, lat, lon]=B1n+P1n+P2n+P3n+P4n+R1n+R6n+Z3n+Z4n+Z5n+Z6n; BPZn[time, depth, lat, lon]=B1n+P1n+P2n+P3n+P4n+Z3n+Z4n+Z5n+Z6n; RRRn[time, depth, lat, lon]=R1n+R6n; BACn[time, depth, lat, lon]=B1n; PHYn[time, depth, lat, lon]=P1n+P2n+P3n+P4n; MZOn[time, depth, lat, lon]=Z3n+Z4n; mZOn[time, depth, lat, lon]=Z5n+Z6n" -v ./ave.TEST02.nc ./N.TEST02.nc

 ncap -O -s "TOTp[time, depth, lat, lon]=B1p+N1p+P1p+P2p+P3p+P4p+R1p+R6p+Z3p+Z4p+Z5p+Z6p; PO4p[time, depth, lat, lon]=N1p; ORGp[time, depth, lat, lon]=B1p+P1p+P2p+P3p+P4p+R1p+R6p+Z3p+Z4p+Z5p+Z6p; BPZp[time, depth, lat, lon]=B1p+P1p+P2p+P3p+P4p+Z3p+Z4p+Z5p+Z6p; RRRp[time, depth, lat, lon]=R1p+R6p; BACp[time, depth, lat, lon]=B1p; PHYp[time, depth, lat, lon]=P1p+P2p+P3p+P4p; MZOp[time, depth, lat, lon]=Z3p+Z4p; mZOp[time, depth, lat, lon]=Z5p+Z6p" -v ./ave.TEST02.nc ./P.TEST02.nc

 ncap -O -s "CHL[time, depth, lat, lon]=P1l+P2l+P3l+P4l" -v ./ave.TEST02.nc ./CHL.TEST02.nc


echo 'FIN'
