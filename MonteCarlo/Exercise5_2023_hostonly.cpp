#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include <chrono>       // get computer chrono (higher resolution than time)
#include "time.h"       // get computer clock time
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <boost/random.hpp>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

class Exercice7
{

private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  double tfin;          // Temps final
  unsigned int nsteps;  // Nombre de pas de temps
  double D;        // Coefficient de diffusion
  double v0;       // Moyenne de la Gaussienne initiale
  double gamma;    // Coefficient de friction
  double vc;       // Vitesse critique pour la friction
  int    N_part;   // Nombre de particules numériques
  int    N_bins;   // Nombre de bins de l'histogramme
  double vlb;      // v_min des bins de l'histogramme
  double vd_D;     // pour la distribution uniforme entre [vg_D et vd_D]
                   // ou double Dirac f = 1/2 (\delta(v-vg_D) + \delta(v-vd_D))
  double vg_D;     // pour la distribution uniforme entre [vg_D et vd_D] 
                   // ou double Dirac f = 1/2 (\delta(v-vg_D) + \delta(v-vd_D))
  string initial_distrib; // type de distribution initiale ('D' pour Dirac)
  double vhb;      // v_max pour les bins de l'histogramme
  double sigma0;   // écart-type de la Gaussienne initiale
  double prefactor; // \sqrt{2 D \Delta t}
  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  void printOut(bool write)
  {
    valarray<double> moments;  
    double var;
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      bins    = binning_fun(v);
      moments = fun_moments(v) ;
      var     = moments[1] - moments[0]*moments[0];
      *outputFile << t << " " << N_part<< " "<<moments[0]<<" "<<var;
      for(int ib = 0; ib < N_bins; ++ib){
        *outputFile<<" "<< bins[ib];
      }
      *outputFile <<endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

protected:
  double t;
  double dt;
  valarray<double> v ;
  valarray<double> bins ;

  double acceleration(double v_p){
    return -gamma*(v_p - vc);
  }

  valarray<double> binning_fun(const valarray<double> &v){
     valarray<double> bins_fun = valarray<double>(N_bins); 
     double h = (vhb - vlb)/N_bins;
     int index;
     for(int ip = 0; ip < N_part; ++ip){
       index = int((v[ip] - vlb)/h);
       if( index >= 0 && index < N_bins ){
        bins_fun[index]=bins_fun[index]+1;
       }
     }
     return bins_fun;
  }


  valarray<double> fun_moments(const valarray<double> &v){
     valarray<double> moments = valarray<double>(2); 
     for(int ip = 0; ip < N_part; ++ip){
        moments[0] += v[ip];
        moments[1] += v[ip] * v[ip];
     }
     moments[0] = moments[0]/N_part; 
     moments[1] = moments[1]/N_part; 
     return moments;
  }

  valarray<double> initialisation(){
     valarray<double> vel = valarray<double>(N_part);
     uint64_t  time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>
       (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
     cout << "time_ns= "<<time_ns<<endl;
     
     //boost::mt19937 rng(time(0));
     boost::mt19937 rng(time_ns); //(time_span);
     cout << " time(0) "<<time(0)<<endl;
     if (initial_distrib == "D"){
	cout << "Delta distribution"<<endl;
        for(int ip = 0; ip < int(N_part/2); ++ip){
          vel[ip] = vg_D;
        }
        for(int ip = int(N_part/2); ip < N_part; ++ip){
          vel[ip] = vd_D;
        }
     }
     else{
	cout << "Gaussian distribution"<<endl;
        boost::normal_distribution<> initial_distribution(v0,sigma0);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > initial_velocities_generator(rng,initial_distribution);
        //boost::uniform_int<> initial_distribution(vg_D,vd_D);
	//boost::variate_generator<boost::mt19937&, boost::uniform_int<> > initial_velocities_generator(rng,initial_distribution);
        for(int ip = 0; ip < N_part; ++ip){
          vel[ip] = initial_velocities_generator();
        }
     }
     return vel;
   }


public:

  Exercice7(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.
    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tfin     = configFile.get<double>("tfin",10);
    nsteps   = configFile.get<unsigned int>("nsteps",100); 
    D        = configFile.get<double>("D",0.1);
    gamma    = configFile.get<double>("gamma",0.08);
    v0       = configFile.get<double>("v0",v0);
    sigma0   = configFile.get<double>("sigma0",1);
    N_part   = configFile.get<double>("N_part",10000);
    N_bins   = configFile.get<double>("N_bins",30);
    vhb      = configFile.get<double>("vhb",5);
    vlb      = configFile.get<double>("vlb",-5);
    vd_D     = configFile.get<double>("vd_D",2);
    vg_D     = configFile.get<double>("vg_D",-2);
    vc       = configFile.get<double>("vc",0);
    sampling = configFile.get<unsigned int>("sampling",2); 
    initial_distrib = configFile.get<string>("initial_distrib","G");

    dt = tfin / nsteps; 

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales

  };

  ~Exercice7()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    bins =valarray<double>(N_bins);
    v =valarray<double>(N_part); 
    prefactor = sqrt(2.*D*dt);
    t = 0.; // initialiser le temps
    last = 0; // initialise le parametre d'ecriture
    boost::mt19937 rng;
    v = initialisation();
    printOut(true); // ecrire premier pas de temps
    boost::normal_distribution<> displace_gauss(0.,1.);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_deplacement(rng,displace_gauss);
    
    for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
    {
      for(int ip = 0; ip < N_part; ++ip){
        v[ip] += dt * acceleration(v[ip]) + prefactor * random_deplacement();
      }
      t = t + dt;
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  };
};


int main(int argc, char* argv[])
{
  Exercice7 engine(argc,argv); // definer la class pour la simulation
  auto start = std::chrono::high_resolution_clock::now();
  engine.run(); // executer la simulation
  auto end  = std::chrono::high_resolution_clock::now();
  auto diff = std::chrono::duration<double>(end-start);
  printf("The duration of the simulation is %g [s] \n", diff );

  cout << "Fin de la simulation." << endl;
  return 0;
}
