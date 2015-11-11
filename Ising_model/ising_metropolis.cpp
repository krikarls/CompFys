#include <armadillo>
#include <iostream>
#include <random>
#include <functional>
#include <ctime>

using namespace arma;

void ising_metropolis(int ii,int D, int N, double T,Mat<double> &lattice,double mag, double eng, Col<double> &E, Col<double> &M,Col<double> &C_v,Col<double> &chi, long int seed);
void rndm_mtrx(int D,Mat<double> &lattice, std::mt19937_64 &generator);
double rnd(std::mt19937_64 &generator);
int rand_pos(int D, std::mt19937_64 &generator);
double lattice_energy(int D,Mat<double> &lattice);

int main()
{
    std::cout << "This program should give good results for lattices up to around 100x100 with current configurations. For bigger lattices small modifications must be made." << std::endl;

    int N = 100000000;              // total cycles
    int D;                         // lattice dimensioin
    int T_points;                   // number of temperature samples
    double T_a;                     // lower temperature
    double T_b;                     // upper temperature

    std::cout << "Read in lattice size:" << std::endl;
    std::cin >> D;
    std::cout << "Give temperature interval limits:" << std::endl;
    std::cin >> T_a >> T_b;
    std::cout << "Number of mesh points for temperature" << std::endl;
    std::cin >> T_points;

    vec T = linspace<vec>(T_a, T_b, T_points);      // array with all temperatures
    vec E = zeros<vec>(T_points);
    vec M = zeros<vec>(T_points);
    vec C_v = zeros<vec>(T_points);
    vec chi = zeros<vec>(T_points);

    Mat<double> lattice(D,D);

    auto seed = std::time(nullptr);
    std::mt19937_64 generator(seed);

    rndm_mtrx(D,lattice, generator);

    ising_metropolis(0,D,N,T(0),lattice,accu(lattice),lattice_energy(D,lattice),E,M,C_v,chi,seed);

    for (int j=0; j < T_points; j++)
        {
            double mag = accu(lattice);
            double eng = lattice_energy(D,lattice);

            ising_metropolis(j,D,N,T(j),lattice,mag,eng,E,M,C_v,chi,seed);

        }

    std::ofstream outFile;
    outFile.open("therm_data.txt", ios::out);
    outFile << D << " " << N << " " << 0 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < T_points; i++) {
               outFile << T(i) << " " << E(i) << " " << M(i) << " " << C_v(i) << " " << chi(i) << endl;
        }
        outFile.close();


} // ***************  END MAIN  ***************

void ising_metropolis(int ii,int D, int N, double T,Mat<double> &lattice,
                        double mag, double eng, Col<double> &E, Col<double> &M,
                        Col<double> &C_v,Col<double> &chi, long int seed)
{

    int start_gather = 1000000;    // cycles before start gathering data
    int measure = 10000;            // how often sample data

    std::mt19937_64 generator(seed);

    // Precomputed exponents
    vec ExP = zeros<vec>(9);
    ExP(4+4) = exp(-4.*2/T);
    ExP(4+2) = exp(-2.*2/T);
    ExP(4+0) = exp(0.*2/T);
    ExP(4-2) = exp(2.*2/T);
    ExP(4-4) = exp(4.*2/T);

    vec avg = zeros<vec>(5);    // storing data samples

    for (int n=0; n < N; n++)
        {
               int i = rand_pos(D,generator);
               int j = rand_pos(D,generator);

               int S = lattice(i,j);
               int dE = lattice(i,j)*(lattice((i+1)%D, j) + lattice(i,(j+1)%D) + lattice((i-1+D)%D,j) + lattice(i,(j-1+D)%D));

               if (ExP(4+dE) > rnd(generator)){ // flip the spin
                   lattice(i,j) = -S;
                   eng += 2*dE;
                   mag -= 2*S;
               }
               if (n > start_gather && (n%measure) == 0){
                           avg(0) += 1;
                           avg(1) += eng;
                           avg(2) += abs(mag);
                           avg(3) += eng*eng;
                           avg(4) += mag*mag;
               }

         }

    E(ii) = (avg(1)/avg(0))/(D*D) ;
    M(ii) = (avg(2)/avg(0))/(D*D) ;
    C_v(ii) = (avg(3)/avg(0)-(avg(1)/avg(0))*(avg(1)/avg(0)))/((D*D)*(T*T));
    chi(ii) = (avg(4)/avg(0)-(avg(2)/avg(0)*(avg(2)/avg(0))))/((D*D)*T);

}

double lattice_energy(int D, Mat<double> &lattice)
{
    double eng = 0;

    for (int i=0; i < D; i++)
        {
        for (int j=0; j < D; j++)
            {
                   eng += -lattice(i,j)*( lattice((i+1)%D,j) + lattice(i,(j+1)%D) + lattice((i-1+D)%D,j) + lattice(i,(j-1+D)%D));
            }
        }
     return eng/2.0; // Each par counted twice
}

int rand_pos(int D, std::mt19937_64 &generator)
{
    std::uniform_int_distribution<int> distribution(0,D-1);
    return distribution(generator);
}

double rnd(std::mt19937_64 &generator)
{
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    return distribution(generator);
}

void rndm_mtrx(int D,Mat<double> &lattice,std::mt19937_64 &generator)
{
    std::uniform_int_distribution<int> distribution(0,1);
    auto RNG = std::bind(distribution,generator);

    for (int ii=0; ii < D; ii++)
    {
       for (int jj=0; jj < D; jj++)
          {
            lattice(ii,jj) = 2*RNG()-1;
          }
    }
}
