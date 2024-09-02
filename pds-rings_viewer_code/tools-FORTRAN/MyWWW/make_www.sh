PATH=/usr/local/bin:/usr/bin gfortran -c *.f
gcc -c *.c
rm -f *.a
ar crs www.a *.o
rm -f *.o
