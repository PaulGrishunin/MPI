#include"SimulationS.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<iomanip>

using namespace std;

SimulationS::SimulationS( long _size, MyMPI *_mmpi ) {
   size = _size;
   initPositions();
   mmpi = _mmpi;
}

void SimulationS::setCumulativeProbability( double **cp ) {
   cumulativeProbability = cp;
}

// polozenie poczatkowe to (0,0)
void SimulationS::initPositions() {
   x = new int[ size ];
   y = new int[ size ];
   
   for ( long i = 0; i < size; i++ ) 
     x[ i ] = y[ i ] = 0;
}

void SimulationS::init() {
   // w wersji sekwencyjnej nie ma tu nic ciekawego
}

long SimulationS::distanceSQ( long i, long j ) {
   long dx = x[ i ] - x[ j ];
   long dy = y[ i ] - y[ j ];
   return dx * dx  + dy * dy;
}

// co sie dzieje z i-ta czastka
int SimulationS::whichOne( long i ) {
   double r = drand48();
   
//   cout << "i " << i << endl;
   for ( int j = 0; j < 5; j++ ) {// stop + 4 kierunki
//      cout << " j " << j << " " << cumulativeProbability[ i ][ j ] << endl;
      if ( r < cumulativeProbability[ i ][ j ] ) return j;
   }
   return 4;
}

void SimulationS::move( long i, int option ) {
   if ( option == 0 ) return; // stoimy w miejscu
   
   switch ( option ) {
     case 1: y[i]++; // w gore
     	break;
     case 2: x[i]++; // w prawo
     	break;
     case 3: y[i]--; // w dol
     	break;
     case 4: x[i]--; // w lewo
     	break;
    }
}

void SimulationS::doIt( int steps  ) {
   for ( int step = 0; step < steps; step++ )
      for ( long i = 0; i < size; i++ )
         move( i, whichOne( i ) );
}

double SimulationS::avgX() {
   double avg = 0;
   for ( long i = 0; i < size; i++ )
     avg += x[i];

   return avg / size;
}

double SimulationS::avgY() {
   double avg = 0;
   for ( long i = 0; i < size; i++ )
     avg += y[i];

   return avg / size;
}

// liczymy średnią liczbę punktów, które rozdziela
// odległość od min do max włącznie.
// tą część kodu można napisać lepiej, ale
// tak będzie łatwiej uzyskać wyższą efektywność.
void SimulationS::calcAvgNumberOfPointsBetween( int min, int max ) {
   long minSQ = min * min;
   long maxSQ = max * max;

   long numberOfPointsBetween = 0;
   long dSQ;

   // dla każdej cząstki
   for ( int i = 0; i < size; i++ ) {
      // przeglądamy jej sąsiadów
      for ( int j = 0; j < size; j++ ) {
         // omijając tu spotkanie cząstki z samą sobą
         if ( i != j ) {
            dSQ = distanceSQ(i, j);
            // zliczamy tych, których odległość zgadza się
            // z zapytaniem

//            cout << "dSQ = " << dSQ << endl;
            if ( ( dSQ >= minSQ ) && ( dSQ <= maxSQ ) ) 
               numberOfPointsBetween++;
         }
      }
   }

   avgNumberOfPointsBetween = (double)numberOfPointsBetween / (double)size;
}

double SimulationS::getAvgNumberOfPointsBetween() {
   return avgNumberOfPointsBetween;
}

ostream & operator<< ( ostream &o, const SimulationS &s ) {
	o << "Points # " << s.size << endl;
	for ( int i = 0; i < s.size; i++ ) {
	  o << setw( 5 ) << i << " " << setprecision(4) << fixed << setw( 7 ) << s.x[i] << " " << right << setw( 7 ) << s.y[ i ] << endl;
	}
	return o;
}
