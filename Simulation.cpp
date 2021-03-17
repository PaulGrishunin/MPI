#include"Simulation.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<iomanip>

using namespace std;

Simulation::Simulation( long _size, MyMPI *_mmpi ) {
   size = _size;
   initPositions();
   mmpi = _mmpi;
}

void Simulation::setCumulativeProbability( double **cp ) {
   cumulativeProbability = cp;                               
}

// polozenie poczatkowe to (0,0)
void Simulation::initPositions() {
   x = new int[ size ];
   y = new int[ size ];
   
   for ( long i = 0; i < size; i++ ) 
     x[ i ] = y[ i ] = 0;
}

void Simulation::init() {                            
 
  int myRank;
    mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    double * buffer = new double[ size * 5 ];
    if ( myRank == 0 ) {
        for ( long i = 0; i < size; i++ )
            for (int j = 0; j < 5; ++j )
                buffer[i*5+j] = cumulativeProbability[i][j];
        mmpi->MPI_Bcast( buffer, 5*size, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }
    else {
        MPI_Status status;
		mmpi->MPI_Bcast( buffer, size * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
        cumulativeProbability = new double * [ size ];
        for ( long i = 0; i < size; i++ ) {
            cumulativeProbability[ i ] = new double [5];
            for (int j = 0; j < 5; ++j )
                cumulativeProbability[i][j] = buffer[i*5+j];
        }
    }
    delete[] buffer;

}
 
 // int myid = mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
 // mmpi->MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  
  // w wersji sekwencyjnej nie ma tu nic ciekawego


long Simulation::distanceSQ( long i, long j ) {
   long dx = x[ i ] - x[ j ];
   long dy = y[ i ] - y[ j ];
   return dx * dx  + dy * dy;    
}

// co sie dzieje z i-ta czastka
int Simulation::whichOne( long i ) {
   double r = drand48();              
   
//   cout << "i " << i << endl;
   for ( int j = 0; j < 5; j++ ) {// stop + 4 kierunki                          
//      cout << " j " << j << " " << cumulativeProbability[ i ][ j ] << endl;
      if ( r < cumulativeProbability[ i ][ j ] ) return j;
   }
   return 4;
}

void Simulation::move( long i, int option ) {
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

void Simulation::doIt( int steps  ) {
    int myRank;
    mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
   int nproc;
   mmpi->MPI_Comm_size(MPI_COMM_WORLD, &nproc);    
   for ( int step = 0; step < steps; step++ ){
      for ( long i = myRank; i < size; i+= nproc ){                
		move( i, whichOne( i ) ); }                             
}
    double sum_x = 0;
	double sum_y = 0;
	
    for ( long i = myRank; i < size; i+= nproc ){
		sum_x += x[i];
		sum_y += y[i];
	}
	
    mmpi->MPI_Reduce( &sum_x, &sumx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mmpi->MPI_Reduce( &sum_y, &sumy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

int flag = 1;
}


double Simulation::avgX() {
 /*  double avg = 0;
   for ( long i = 0; i < size; i++ )
     avg += x[i]; */

   return sumx / size;
}

double Simulation::avgY() {
 /*  double avg = 0;
   for ( long i = 0; i < size; i++ )
     avg += y[i]; */

   return sumy / size;

int flag = 0;
}

// liczymy średnią liczbę punktów, które rozdziela
// odległość od min do max włącznie.
// tą część kodu można napisać lepiej, ale
// tak będzie łatwiej uzyskać wyższą efektywność.

void Simulation::calcAvgNumberOfPointsBetween( int min, int max ) {
   long minSQ = min * min;
   long maxSQ = max * max;

   long numberOfPointsBetween = 0;
   long dSQ;
   

   int myRank;
   mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
   int nproc;
   mmpi->MPI_Comm_size( MPI_COMM_WORLD, &nproc );

    if ( flag == 1 ){
    double * tableX = new double[ size ];
    double * tableY = new double[ size ];
    
	
    for ( int p = 0; p < nproc; p++ ){
	    if ( myRank == p ) {
            for ( long i = myRank; i < size; i+= nproc ){
				tableX[i] = x[i];
				tableY[i] = y[i];
			}
		}
		mmpi->MPI_Bcast( tableX, size, MPI_DOUBLE, p, MPI_COMM_WORLD );
		mmpi->MPI_Bcast( tableY, size, MPI_DOUBLE, p, MPI_COMM_WORLD );
		
		if ( myRank != p ) {
			 for ( long i = p; i < size; i+= nproc ){
				x[i] = tableX[i];
				y[i] = tableY[i];
				
		}
	}
    }
}

    long calc = 0;
   // dla każdej cząstki
   for ( int i = nproc; i < size; i+=nproc ) {
      // przeglądamy jej sąsiadów
      for ( int j = 0; j < size; j++ ) {
         // omijając tu spotkanie cząstki z samą sobą
         if ( i != j ) {
            dSQ = distanceSQ(i, j);
            // zliczamy tych, których odległość zgadza się
            // z zapytaniem
            if ( ( dSQ >= minSQ ) && ( dSQ <= maxSQ ) ) 
               numberOfPointsBetween++;
         }
      }
   } 
		 
   mmpi->MPI_Reduce( &calc, &numberOfPointsBetween, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   
   if ( myRank == 0 ){
   avgNumberOfPointsBetween = (double)calc / (double)size;
   }
}

double Simulation::getAvgNumberOfPointsBetween() {
   return avgNumberOfPointsBetween;
}

ostream & operator<< ( ostream &o, const Simulation &s ) {
	o << "Points # " << s.size << endl;
	for ( int i = 0; i < s.size; i++ ) {
	  o << setw( 5 ) << i << " " << setprecision(4) << fixed << setw( 7 ) << s.x[i] << " " << right << setw( 7 ) << s.y[ i ] << endl;
	}
	return o;
}
