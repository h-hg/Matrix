Demo.exe:Fraction.o Demo.o
	g++ -o Demo.exe Fraction.o Demo.o
fraction.o:Fraction.cpp Fraction.h
	g++ -c Fraction.cpp
Demo.o:Demo.cpp Matrix.hpp
	g++ -c Demo.cpp