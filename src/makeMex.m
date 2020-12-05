
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if isOctave
    # This seems to work on Octave
    mex genUDCPDMex.c udcpd.c misc.c
    mex genVDCPDMex.c udcpd.c vdcpd.c misc.c
else
    mex -DMEX_COMPILE_FLAG -O CFLAGS="\$CFLAGS -std=c99 -pedantic" genUDCPDMex.c udcpd.c misc.c
    mex -DMEX_COMPILE_FLAG -O CFLAGS="\$CFLAGS -std=c99 -pedantic" genVDCPDMex.c udcpd.c vdcpd.c misc.c
end
