FC = gfortran
CC = gcc
driver.exe: driver.f95 adj.o tlm.o PUSHPOPDiff.o adStack.o adBuffer.o
	gfortran -o driver.exe -freal-4-real-16 adj.o tlm.o PUSHPOPDiff.o adStack.o adBuffer.o -freal-4-real-16 driver.f95
adj.o: adj.f95
	gfortran -c adj.f95 -freal-4-real-16
tlm.o: tlm.f95
	gfortran -c tlm.f95 -freal-4-real-16
PUSHPOPDiff.o: PUSHPOPDiff.f
	gfortran -c PUSHPOPDiff.f -freal-4-real-16
adStack.o: adStack.c
	gcc -c adStack.c
adBuffer.o: adBuffer.f
	gfortran -c adBuffer.f -fallow-argument-mismatch -freal-4-real-16
clean:
	rm *.o
