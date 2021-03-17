#ifndef SIMULATION_H

#define SIMULATION_H

#include<iostream>
#include"MyMPI.h"

using namespace std;

class Simulation {
	private:
		MyMPI *mmpi;
		double avgNumberOfPointsBetween;
		double sumx;
		double sumy;
		int flag = 0;          // флаг показывает был ли выполнен метод doIt (1) или avg (0) 
		
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

		Simulation( long _size, MyMPI *_mmpi );

	friend  ostream & operator<<( ostream &o, const Simulation &s );
};

#endif
