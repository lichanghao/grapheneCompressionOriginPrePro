 ifort -c -w95 -O3 headers.f90 
 ifort -c -w95 -O3 *.f90 
 ifort -O3 -o PrePro *.o
 rm -rf *.o

