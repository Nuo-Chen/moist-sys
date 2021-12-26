FC = gfortran
CC = gcc
driver.exe: driver.f ode5123_fwd_b.o  ode5123_fwd_d.o PUSHPOPDiff.o adStack.o adBuffer.o
	gfortran -o driver.exe -freal-4-real-16 ode5123_fwd_b.o ode5123_fwd_d.o PUSHPOPDiff.o adStack.o adBuffer.o -freal-4-real-16 driver.f 
ode5123_fwd_b.o: ode5123_fwd_b.f
	gfortran -c ode5123_fwd_b.f -freal-4-real-16
ode5123_fwd_d.o: ode5123_fwd_d.f
	gfortran -c ode5123_fwd_d.f -freal-4-real-16
PUSHPOPDiff.o: PUSHPOPDiff.f
	gfortran -c PUSHPOPDiff.f -freal-4-real-16
adStack.o: adStack.c
	gcc -c adStack.c
adBuffer.o: adBuffer.f
	gfortran -c adBuffer.f -Wno-argument-mismatch -freal-4-real-16
clean:
	rm *.o
