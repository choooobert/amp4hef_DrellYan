#include <iostream>
#include <complex>
#include "amp4hef.hpp"
#include <cmath>
#define pi 3.1415926535897932384626433832795
#include <string>
#include <exception>
#include <fstream>
#include <cstdio>
#include <vector>


using namespace std;


void read_event( ifstream &file ,bool &exitLoop ,double &eventWeight
                ,vector<double> &momenta ,const int &Nfinst
                ,double &instWeight ,double &psWeight ,double &partonLumi
                ,double &alphaStrong );

double sHat(vector <double> momenta);

int main(int argc, char * argv[]){

      bool exitLoop = false;
      int Noffshell,  Nfinst, Ntotal ,itmp ,id1;
      double Ecm, kTsq[2];
      double eventWeight, instWeight, psWeight, cnstWeight, totalWeight;
      double partonLumi, alphaStrong, flux ,ampSquared;
      double sumW0 = 0, sumW1 = 0, sumW2 = 0;
      vector <int> process;
      string eventFile = argv[1];
      ifstream file;

// Read the center-of-mass energy, the number of off-shell partons,
// the number of final-state partons, and the process.
      try{
            file.open(eventFile);
            string dane;
            if (!file.is_open())
                  throw string("error opening file\n");

            do{
                  getline(file, dane);
            } while (dane != "BEGINOFFILE");

            file >> Ecm;
            file >> Noffshell;
            file >> Nfinst;
            getline(file, dane);
            for (int i = 0; i < Nfinst+2; ++i){
                  file >> itmp;
                  process.push_back(itmp);
            }
            getline(file, dane);
            Ntotal = Nfinst + 2;
      }
      catch (string state){
            cout << state;
            return 0;
      }

// Put the processes, and get id.
      put_process( &id1 ,&Ntotal ,&Noffshell ,&process.at(0) );

// Calculate the overall constant to the event weights.
      cnstWeight = 1
//     Conversion factor from GeV to nanobarn.
                 * 389379.66 
//     The phase space weight "psWeight" is missing the following factor,
//     according to the RAMBO convention.
                 * pow(2*pi, 4-3*Nfinst)
//     Average over initial-state spins. The weight factor "partonLumi" includes a
//     factor 2 for each off-shell gluon to have a uniform definition of cnstWeight.
                 / (2 * 2)
//     Conversion factor from alphaStrong to g_QCD^2
                 * pow(4*pi, Nfinst);
//     Average over initial-state colors and symmetry factor are included in ampSquared.

      vector <double> momenta;

      double directions[2 * 4];

      while (true){
      
// Besides the external momenta, read_event reads the following from the event file:
//  instWeight: the weight from the generation of initial-state variables
//    psWeight: the weight from the generation of final-state phase space
//  partonLumi: product of the pdfs for the initial-state partons
// alphaStrong: the value of the strong coupling
// eventWeight: the correct value of the weight determined during the
//                creation of the event file
            read_event( file ,exitLoop ,eventWeight ,momenta ,Nfinst
                       ,instWeight ,psWeight ,partonLumi ,alphaStrong );
            if (exitLoop)
                  break;

// The event file includes events that fell outside the cuts.
// These count as events where the integrand vanishes.
            ++sumW0;

            if (eventWeight <= 0E0)
                  continue;

// Determine the flux factor.
            kTsq[0] = momenta[1]*momenta[1] + momenta[2]*momenta[2];
            kTsq[1] = momenta[5]*momenta[5] + momenta[6]*momenta[6];
            flux = pow(sHat(momenta) + kTsq[0] + kTsq[1], 2) - 4*kTsq[0]*kTsq[1];
            flux = 2*sqrt( flux );

// Construct directions from the energy and the z-component of the initial-state momenta.
            directions[0]=momenta[0];  directions[4]=momenta[4];
            directions[1]=0;           directions[5]=0;
            directions[2]=0;           directions[6]=0;
            directions[3]=momenta[3];  directions[7]=momenta[7];

// Evaluate the matrix element. Call matrix_element_b which includes factors to
// average over initial-state colors, and the final-state symmetry factor.
            put_momenta( &id1, &momenta.at(0), directions );
            matrix_element_b( &id1 ,&ampSquared );

// Determine the total weight of the event.
// Choose ampSquared[1] which is averaged over initial-state colors and
// includes the final-state symmetry factor.
            totalWeight = cnstWeight / flux * instWeight * psWeight * partonLumi
                        * ampSquared * pow(alphaStrong,Nfinst);

// Gather statistics.
            sumW1 = sumW1 + totalWeight;
            sumW2 = sumW2 + totalWeight*totalWeight;

// Compare calculated weight with the number from the file.
            printf( "( re-calculated weight )/( weight from file ): %20.19f \n"
                   , totalWeight / eventWeight );
            momenta.clear();
      }


// Finalize statistics.
      cout << "\ncross section (in nb): "
           << sumW1 / sumW0 << endl;
      cout << "error estimate (in %):" 
           << 100 * sqrt((sumW2*sumW0 / pow(sumW1, 2) - 1) / (sumW0 - 1)) << endl;

}


void read_event( ifstream &file ,bool &exitLoop ,double &eventWeight
                ,vector <double> &momenta ,const int &Nfinst
                ,double &instWeight ,double &psWeight ,double &partonLumi
                ,double &alphaStrong ){
      string tmp;
      double tmp2;
      getline(file, tmp);
      
      if (tmp == "ENDOFFILE") {
            
            exitLoop = true;
            return;
      }

      else if (tmp.find("EVENT") != std::string::npos) {
            
            file >> eventWeight;
            
            if (eventWeight < 10E-15) {
                  getline(file, tmp);
                  return;
            }

            for (int i = 0; i < 4 * 2 + 4 * Nfinst; ++i){
                  file >> tmp2;
                  momenta.push_back(tmp2);
            }
            file >> instWeight;
            file >> psWeight;
            file >> partonLumi;
            file >> alphaStrong;
            getline(file, tmp);
      }


}

double sHat(vector <double> momenta){
      return( (momenta[0]+momenta[4]+momenta[3]+momenta[7])
             *(momenta[0]+momenta[4]-momenta[3]-momenta[7])
             -(momenta[1] + momenta[5])*(momenta[1] + momenta[5])
             -(momenta[2] + momenta[6])*(momenta[2] + momenta[6]) );
}
