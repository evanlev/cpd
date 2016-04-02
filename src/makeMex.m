mex -DMEX_COMPILE_FLAG -O CFLAGS="\$CFLAGS -std=c99 -pedantic" genUDCPDMex.c udcpd.c misc.c
mex -DMEX_COMPILE_FLAG -O CFLAGS="\$CFLAGS -std=c99 -pedantic" genVDCPDMex.c udcpd.c vdcpd.c misc.c

