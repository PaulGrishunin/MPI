#ifndef SIMULATIONS_H

#define SIMULATIONS_H

#include<iostream>
#include"MyMPI.h"

using namespace std;

class SimulationS {
	private:
		MyMPI *mmpi;
		double avgNumberOfPointsBetween;
	protected:
		int *x, *y; // polozenia punktow
		long size;   // liczba punktow
		double **cumulativeProbability;
		
		long distanceSQ( long i, long j ); // zwraca kwadrat odleglosci pomiedzy i oraz j
		int whichOne( long i ); // zwraca wybrana opcje 0 - stop, 1-4 kierunki ruchu
		void move( long i, int option ); // przesuwa punkt i wg. wskazanej opcji  
		void calcProbablity( double *c ); // przelicza wspolczynniki na prawdopodobienstwa
        void initPositions();

	public:
		void setCumulativeProbability( double **cp ); 
		void init();

		void doIt( int steps ); // ta metoda zleca wykonanie obliczen
		double avgX();
		double avgY();
		void calcAvgNumberOfPointsBetween( int min, int max ); 
		double getAvgNumberOfPointsBetween(); // pobranie wyniku zleconego metodą powyzej

		SimulationS( long _size, MyMPI *_mmpi );

	friend  ostream & operator<<( ostream &o, const SimulationS &s );
};

#endif
