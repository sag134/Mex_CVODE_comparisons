#!/bin/bash
gcc -Wall -O3 MichMent_gcc_looped.c -o MichMent_gcc_looped.exe /home/sag134/Sundials/sundials-2.6.2/instdir/lib/libsundials_cvode.a /home/sag134/Sundials/sundials-2.6.2/instdir/lib/libsundials_nvecserial.a -lm 
./MichMent_gcc_looped.exe > out

