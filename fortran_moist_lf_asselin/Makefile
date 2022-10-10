FC = gfortran 
CC = gcc
CFLAGS = -target x86_64-apple-darwin20.3.0
driver.exe: driver.f ode_lf_asselin_adj.o ode_lf_asselin_tlm.o PUSHPOPDiff.o adStack.o adBuffer.o
	gfortran -o driver.exe -freal-4-real-16 ode_lf_asselin_adj.o ode_lf_asselin_tlm.o PUSHPOPDiff.o adStack.o adBuffer.o -freal-4-real-16 driver.f
ode_lf_asselin_adj.o: ode_lf_asselin_adj.f
	gfortran -c ode_lf_asselin_adj.f -freal-4-real-16
ode_lf_asselin_tlm.o: ode_lf_asselin_tlm.f
	gfortran -c ode_lf_asselin_tlm.f -freal-4-real-16
PUSHPOPDiff.o: PUSHPOPDiff.f
	gfortran -c PUSHPOPDiff.f -freal-4-real-16
adStack.o: adStack.c
	gcc -c adStack.c
adBuffer.o: adBuffer.f
	gfortran -c adBuffer.f -fallow-argument-mismatch -freal-4-real-16
clean:
	rm *.o