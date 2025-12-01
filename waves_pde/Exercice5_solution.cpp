#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.h"
#include <algorithm>

using namespace std;



void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, \
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
      if (bc_l == "fixed"){
        fnext[0] = 0; // TODO : Completer la condition au bord gauche fixe
      }else if(bc_l == "free"){
        fnext[0] = fnext[1]; // TODO : Completer la condition au bord gauche libre
      }else if (bc_l =="exit"){
        fnext[0] = fnow[0] + (fnow[1]-fnow[0])*sqrt(beta2[0]); // TODO : Completer la condition au bord gauche "sortie de l'onde"
      }else{
        cerr << "Merci de choisir an condition aux bord valid pour la gauche" << endl;
      }
	      
      if (bc_r == "fixed"){
        fnext[N-1] = 0; // TODO : Completer la condition au bord droit fixe
      }else if(bc_r == "free"){
        fnext[N-1] = fnext[N-2]; // TODO : Completer la condition au bord droit libre
      }else if (bc_r =="exit"){
        fnext[N-1] = fnow[N-1] + (fnow[N-2] - fnow[N-1]) * sqrt(beta2[N-1]); // TODO : Completer la condition au bord droit "sortie de l'onde"
      }else{
        cerr << "Merci de choisir an condition aux bord valid pour la droit" << endl;
      }
}

double finit(double x, double xL, double n_init, double xR, double A, double x1, double x2, string initialization)
{
  double finit_(0.);
  const double PI = 3.1415926535897932384626433832795028841971e0;

if(initialization=="mode"){
// TODO: initialiser un mode propre
  finit_ = sin((n_init + 0.5) * PI * (x - xL) /(xR-xL));
}
else{
  finit_ = A*0.5*(1-cos(2.0 * PI * (x - x1)/(x2 -x1)))*(x1<x && x<x2) ;
}
  return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0), N;

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Simulation parameters :
  int nx                    = configFile.get<int>("nx",64); // nb intervalles
  int n_stride              = configFile.get<int>("n_stride",10);

  double tfin                  = configFile.get<double>("tfin",10);
  double CFL                   = configFile.get<double>("CFL",1);
  double nsteps                = configFile.get<double>("nsteps",3);
  double A                     = configFile.get<double>("A",1);
  double n_init                = configFile.get<double>("n_init",3);
  double hL                    = configFile.get<double>("hL",7000);
  double hR                    = configFile.get<double>("hR",200);
  double hC                    = configFile.get<double>("hC",35);
  double h00                   = configFile.get<double>("h00",3);
  double x1                    = configFile.get<double>("x1",2.0);
  double x2                    = configFile.get<double>("x2",6.0);
  double xa                    = configFile.get<double>("xa",3*pow(10,5));
  double xb                    = configFile.get<double>("xb",7*pow(10,5));
  double xc                    = configFile.get<double>("xc",7.2*pow(10,5));
  double xd                    = configFile.get<double>("xd",8.5*pow(10,5));
  double xL                    = configFile.get<double>("xL",0);
  double xR                    = configFile.get<double>("xR",10.0);
  bool impose_nsteps           = configFile.get<bool>("impose_nsteps",false);
  bool v_uniform               = configFile.get<bool>("v_uniform",true);
  bool usecub                  = configFile.get<bool>("usecub",true);
  bool bc_cuda                 = configFile.get<bool>("bc_cuda",true);
  bool swap_bool               = configFile.get<bool>("swap_bool",true);
  bool thrust_device           = configFile.get<bool>("thrust_device",true);
  bool ecrire_f                = configFile.get<bool>("ecrire_f",true);
  string equation_type         = configFile.get<string>("equation_type","Eq1");
  string bc_l              = configFile.get<string>("cb_gauche","free");
  string bc_r              = configFile.get<string>("cb_droite","free");
  string initialization        = configFile.get<string>("initialization","mode");
  string initial_state         = configFile.get<string>("initial_state","static");
  string output                = configFile.get<string>("output","outputout");

  N = nx+1;                                // pts of the mesh

  vector<double> h0(N) ;
  vector<double> vel2(N) ;
  vector<double> x(N) ;
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = (xR - xL) / (N-1);
  

  for(int i(0); i<N; ++i){ 
     x[i] = xL + i * dx ;
     h0[i] = 0.0;
     if(v_uniform){
           h0[i]  = h00;
     } 
     else {
           h0[i]  = hL * (xL<=x[i] && x[i]<=xa) + \
		   (0.5*(hL + hC) + 0.5*(hL - hC)*cos(PI *(x[i]-xa)/(xb - xa))) * (xa<x[i] && x[i]<xb) +\
          hC * (xb<=x[i] && x[i]<=xc) + \
		   (0.5*(hR + hC) - 0.5*(hR - hC)*cos(PI *(x[i]-xc)/(xd - xc))) * (xc<x[i] && x[i]<xd) +\
		     hR * (xd<=x[i] && x[i]<=xR);


	   if(i==0) cout << "v is not uniform"<<endl;
     }
     vel2[i]  = g* h0[i];
  }

  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());
  dt = CFL * dx / sqrt(*max_vel2);

  if(impose_nsteps){
    dt  = tfin/nsteps;
    CFL = dt/dx * sqrt(*max_vel2);
  }


  ofstream fichier_x((output + "_x").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f").c_str());
  fichier_f.precision(15);

  // Initialisation des tableaux du schema numerique :

  //TODO initialize f and beta
  for(int i(0); i<N; ++i)
  {
    fpast[i] = 0.;
    fnow[i]  = 0.;
    beta2[i] = vel2[i] * dt * dt / (dx * dx);

    fnow[i]  = finit(x[i],  xL,  n_init,  xR, A, x1, x2, initialization);

    if(initial_state =="static"){
      fpast[i] = fnow[i]; // system is at rest for t<=0, 
    }
    else if(initial_state =="right"){ 
      fpast[i] = finit(x[i] +  sqrt(vel2[i]) * dt, xL,  n_init,  xR, A, x1, x2, initialization);
    }
    else if(initial_state =="left"){
      fpast[i] = finit(x[i] -   sqrt(vel2[i]) * dt, xL,  n_init,  xR, A, x1, x2, initialization);
    }
  }


  cout<<"beta2[0] is "<<beta2[0]<<endl;
  cout<<"dt is "<< dt <<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
     }
    ++stride;

    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
      fnext[i] = 2.*(1.-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+fnow[i-1]); // TODO : Compléter le schéma A
      if(equation_type == "Eq1"){
        fnext[i] += .25*(beta2[i+1]-beta2[i-1])*(fnow[i+1]-fnow[i-1]); // TODO : Compléter le schéma B
      }else if(equation_type == "Eq6"){  // part facultatif
        fnext[i] += .5*(beta2[i+1]-beta2[i-1])*(fnow[i+1]-fnow[i-1]) + (beta2[i+1]-2.*beta2[i]+beta2[i-1])*fnow[i]; // TODO : Compléter le schéma C
      }

    }

    // TODO add boundary conditions
    boundary_condition(fnext, fnow, A, t, dt, beta2, bc_l, bc_r, N);

    // Mise a jour :
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();

  return 0;
}
