gfortran -c carray.f90
g++ -c test_tensor.cc
g++ -o test carray.o test_tensor.o -lgfortran
