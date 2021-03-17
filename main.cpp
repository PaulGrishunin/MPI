#include<iostream>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include"SimulationS.h"
#include"Simulation.h"

using namespace std;

void calcCummulative( double *probability ) {
    double sum = 0.0;
    for ( int i = 0; i < 5; i++ )
        sum += probability[ i ];
    
    for ( int i = 0; i < 5; i++ )
        probability[ i ] /= sum;

    probability[ 1 ] += probability[ 0 ];
    probability[ 2 ] += probability[ 1 ];
    probability[ 3 ] += probability[ 2 ];
    probability[ 4 ] += probability[ 3 ];
}

void copy( double *src, double *dst ) {
    for ( int i = 0; i < 5; i++ ) {
        dst[ i ] = src[ i ];
    }
}

void fill( double *pattern, double **destination, long idxFirst, long idxLast ) {
    for ( long i = idxFirst; i < idxLast; i++ ) {
        copy( pattern, destination[ i ] );
    }
}

double **createPattern( int main, int margin ) {
        double pattern1[ 5 ] = { 0, 1.0, 0.0, 0.0, 0.0 }; 
        double pattern2[ 5 ] = { 0, 0.0, 1.0, 0.0, 0.0 }; 
        double pattern3[ 5 ] = { 0, 0.0, 0.0, 1.0, 0.0 }; 
        double pattern4[ 5 ] = { 0, 0.0, 0.0, 0.0, 1.0 }; 

        calcCummulative( pattern1 );
        calcCummulative( pattern2 );
        calcCummulative( pattern3 );
        calcCummulative( pattern4 );

    int particles = main + margin;

// alokacja pamięci
    double ** probability4particles = new double * [ particles ]; 
    for ( long i = 0; i < particles; i++ ) {
        probability4particles[ i ] = new double [5];
    }

    int chunk = main / 4;

// przygotowanie prawdopodobienstw ruchu cząstek
    fill( pattern1, probability4particles, 0, chunk );
    fill( pattern2, probability4particles, chunk, 2*chunk );
    fill( pattern3, probability4particles, 2*chunk, 3*chunk );
    fill( pattern4, probability4particles, 3*chunk, 4*chunk );
    fill( pattern2, probability4particles, 4*chunk, particles );
/*    
    for ( long i = 0; i < particles; i++ ) {
        cout << "P(" << i << ") = [" ;
        for ( int j = 0; j < 5; j++ )
          cout << probability4particles[i][j] << "," ;
        cout << "]" << endl;
    }
*/
    return probability4particles;
}

void serialSimulation( int main, int margin, int repetitions, 
                       int stepsPerRepetision, 
                       int internalRepetitions, MyMPI *mmpi ) {

    int particles = main + margin;
    int steps = stepsPerRepetision / internalRepetitions;
    int totalSteps = 0;

    ofstream report;

    SimulationS *s = new SimulationS( particles , mmpi );

	int myRank;
	mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

    if ( myRank == 0 ) {
        report.open( "serial.txt" );
        s->setCumulativeProbability( createPattern(main,margin) );
    }

    double v1, v2;
    double dT_sum, nT_sum, total_sum;
    double dT_start, nT_start;

    dT_sum = nT_sum = total_sum = 0.0;

    s->init();
    for ( int i = 0; i < repetitions; i++ ) {

        if ( myRank == 0 ) {
            dT_start = mmpi->MPI_Wtime();
        }

        for ( int j = 0; j < internalRepetitions; j++) {
            s->doIt( steps );      
            totalSteps += steps;
            if ( myRank == 0 ) {
                cout << setw(6) << totalSteps ;
                cout << " <x> = " << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgX() 
                     << " <y> = " << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgY()
                     << endl; 
                report << setw(6) << totalSteps ;
                report << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgX() 
                       << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgY()
                     << endl; 
            }
        }

        if ( myRank == 0 ) {
            dT_sum += mmpi->MPI_Wtime() - dT_start;
            nT_start = mmpi->MPI_Wtime();
        }

        s->calcAvgNumberOfPointsBetween(0,32000);
        
        if ( myRank == 0 ) {
            v1 = s->getAvgNumberOfPointsBetween();
        }
        
        s->calcAvgNumberOfPointsBetween(32001,46000);

        if ( myRank == 0 ) {
            nT_sum += mmpi->MPI_Wtime() - nT_start;
            total_sum += mmpi->MPI_Wtime() - dT_start;

            v2 = s->getAvgNumberOfPointsBetween();
            cout << setw(6) << totalSteps ;
            cout << " v1 = "  << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v1
                 << " v2 = "  << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v2 << endl;
            report << setw(6) << totalSteps ;
            report  << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v1
                    << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v2 << endl;
        }
    }

    if ( myRank == 0 ) {
        cout << setw(7) << setprecision(3);
        cout << "doIt() avg time .......................  " << dT_sum / repetitions << endl;
        cout << "calcAvgNumberOfPointsBetween() avg time  " << nT_sum / repetitions << endl;
        cout << "Total avg time ........................  " << total_sum / repetitions << endl;

        report << dT_sum / repetitions << endl;
        report << nT_sum / repetitions << endl;
        report << total_sum / repetitions << endl;
        report.close();
    }
}

bool compare( double expected, double actual, double prec, string txt ) {
    double v = fabs( actual - expected );
    if ( v > prec ) {
        cout << "BŁĄD: " << txt << " expected=" << expected << " actual=" << actual << endl;
        return false;
    }
    return true;
}

bool parallelSimulation( int main, int margin, int repetitions, 
                       int stepsPerRepetision, 
                       int internalRepetitions, MyMPI *mmpi, double efficiency ) {

    int particles = main + margin;
    int steps = stepsPerRepetision / internalRepetitions;
    int totalSteps = 0;
    bool ok = true;

    ifstream report;

    Simulation *s = new Simulation( particles , mmpi );

	int myRank;
	mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

    if ( myRank == 0 ) {
        report.open( "serial.txt" );
        s->setCumulativeProbability( createPattern(main,margin) );
    }

    double v1, v2;
    double dT_sum, nT_sum, total_sum;
    double dT_start, nT_start;

    dT_sum = nT_sum = total_sum = 0.0;

    s->init();

    int totalStepsS;
    double axS, ayS;
    double v1S, v2S;

    for ( int i = 0; i < repetitions; i++ ) {

        if ( myRank == 0 ) {
            dT_start = mmpi->MPI_Wtime();
        }

        for ( int j = 0; j < internalRepetitions; j++) {
            s->doIt( steps );      
            if ( myRank == 0 ) {
                totalSteps += steps;
                cout << setw(6) << totalSteps ;
                cout << " <x> = " << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgX() 
                     << " <y> = " << setiosflags (ios::showbase) << setw(7) << setprecision(3) << s->avgY()
                     << endl; 

                report >> totalStepsS;
                report >> axS;
                report >> ayS;

                ok &= compare( totalStepsS, totalSteps, 0.1, "Liczba kroków symulacji");
                ok &= compare( axS, s->avgX(), 0.01, "Srednia X");
                ok &= compare( ayS, s->avgY(), 0.01, "Srednia Y");
            }
        }

        if ( myRank == 0 ) {
            dT_sum += mmpi->MPI_Wtime() - dT_start;
            nT_start = mmpi->MPI_Wtime();
        }

        s->calcAvgNumberOfPointsBetween(0,32000);
        
        if ( myRank == 0 ) {
            v1 = s->getAvgNumberOfPointsBetween();
        }
        
        s->calcAvgNumberOfPointsBetween(32001,46000);

        if ( myRank == 0 ) {
            nT_sum += mmpi->MPI_Wtime() - nT_start;
            total_sum += mmpi->MPI_Wtime() - dT_start;

            v2 = s->getAvgNumberOfPointsBetween();
            cout << setw(6) << totalSteps ;
            cout << " v1 = "  << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v1
                 << " v2 = "  << setiosflags (ios::showbase | ios::showpoint | ios::fixed) << setw(8) << setprecision(1) << v2 << endl;

            report >> totalStepsS;
            report >> v1S;
            report >> v2S;

            ok &= compare( totalStepsS, totalSteps, 0.1, "Liczba kroków symulacji");
            ok &= compare( axS, s->avgX(), 0.1, "v1");
            ok &= compare( ayS, s->avgY(), 0.1, "v2");
        }
    }

    if ( myRank == 0 ) {
        cout << setw(7) << setprecision(3);
        double dTS, cAS, taS;
        double dTP, cAP, taP;

        dTP = dT_sum / repetitions;
        cAP = nT_sum / repetitions;
        taP = total_sum / repetitions;

        report >> dTS;
        report >> cAS;
        report >> taS;

        double dTprzysp, cAprzysp, taprzysp;
        dTprzysp = dTS / dTP;
        cAprzysp = cAS / cAP;
        taprzysp = taS / taP;

      	int size;
        mmpi->MPI_Comm_size( MPI_COMM_WORLD, &size );
        double dTE, cAE, taE;

        dTE = 100.0*dTprzysp / size;
        cAE = 100.0*cAprzysp / size;
        taE = 100.0*taprzysp / size;

        cout << "doIt() avg time .......................  " << dTP 
            << " ( " << dTprzysp << "x / " << dTE << "100% )" << endl;
        cout << "calcAvgNumberOfPointsBetween() avg time  " << cAP 
            << " ( " << cAprzysp << "x / " << cAE << "100% )" << endl;
        cout << "Total avg time ........................  " << taP 
            << " ( " << taprzysp << "x / " << taE << "100% )" << endl;

        if ( dTE < efficiency ) {
            ok = false;
            cout << "Uzyskano za małą efektywność pracy programu w metodzie doIT, oczekiwano " << efficiency << endl;
        }
        if ( cAE < efficiency ) {
            ok = false;
            cout << "Uzyskano za małą efektywność pracy programu w metodzie calcAvgNumberOfPointsBetween, oczekiwano " << efficiency << endl;
        }

        report.close();
    }

    return ok;
}


int main( int ac, char **av ) {
    long main = 60000   ;
    long margin = 7;
    int repetitions = 10;
    int internalRepetitions = 3;
    int stepsPerRepetision = 9000;

	MyMPI *mmpi = new MyMPI();
	mmpi->MPI_Init( &ac, &av );
  	int rank;
    mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    int canExecuteTest = 1;

    if ( rank == 0 ) {
        ifstream report( "serial.txt");
        if ( ! report.good() ) {
            report.close();
            cout << "Brak pliku z sekwencyjnym wykonaniem zadania!" << endl;
          	int size;
            mmpi->MPI_Comm_size( MPI_COMM_WORLD, &size );
            if ( size > 1 ) {
                cout << "Proszę uruchomić program z jednym procesem!" << endl;
                canExecuteTest = 0;
            } else {
                serialSimulation( main, margin, repetitions, 
                          stepsPerRepetision, internalRepetitions, mmpi );
            }
        } else {
            cout << "Znaleziono plik z zapisem wykonania sekwenyjnego" << endl;
        }
    }

    MPI_Bcast( &canExecuteTest, 1, MPI_INT, 0, MPI_COMM_WORLD );

    if ( canExecuteTest ) {
        cout << "TEST" << endl;
        bool result = parallelSimulation( main, margin, repetitions, 
                            stepsPerRepetision, internalRepetitions, mmpi, 70.0 );
        if ( rank == 0 ) {
            if ( result ) {
                cout << "######  TEST ZALICZONY" << endl;
            } else {
                cout << "######  TEST ZAKONCZONY BLEDEM" << endl;
            }
        }
    }

	mmpi->MPI_Finalize();
    return 0;
}

