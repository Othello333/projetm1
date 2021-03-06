/*
!Ce code C++ modelise le mouvement d'une planete autour du soleil. GM=1, M>>m

Compilation:
c++ -Wall -Wextra equadiff.cpp -o equadiffcpp.x

Temps pris par la fonction principale: 297211 microsecondes, soit 0.297211 seconde

Gnuplot:
plot 'equadiffrk4cpp.res' using 2:3
set xlabel 'x'
set ylabel 'y'
replot

_________________________

X(1)=x(t)
X(2)=y(t)
X(3)=vx(t)
X(4)=vy(t)

Xderiv(1)=vx(t)
Xderiv(2)=vy(t)
Xderiv(3)=-x/(r^3)
Xderiv(4)=-y/(r^3)

conditions initiales annoncees dans l'enonce // utilisees ici
x0=0 // x0=1
y0=1 // y0=0
vx0=-0.5 // vx0=0
vy0=0 // vy0=-0.5
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include <chrono> //pour mesurer le temps d'execution

using namespace std;
using namespace std::chrono;

void deriv(double t, double X[], double Xderiv[], int n){
double radius;
radius=sqrt(pow(X[0],2)+pow(X[1],2));
Xderiv[0]=X[2];
Xderiv[1]=X[3];
Xderiv[2]=-X[0]/pow(radius,3);
Xderiv[3]=-X[1]/pow(radius,3);
}



void rk4(double t, double X[], double dt, int n, void deriv(double, double[], double[], int)){

int i;
double ddt;
double Xp[n], k1[n], k2[n], k3[n], k4[n];
ddt=0.5*dt;

deriv(t,X,k1,n);

for(i=0;i<n;i++){
Xp[i]=X[i]+ddt*k1[i];
deriv(t+ddt,Xp,k2,n);
}

for(i=0;i<n;i++){
Xp[i]=X[i]+ddt*k2[i];
deriv(t+ddt,Xp,k3,n);
}

for(i=0;i<n;i++){
Xp[i]=X[i]+dt*k3[i];
deriv(t+dt,Xp,k4,n);
}

for(i=0;i<n;i++){
X[i] = X[i] + dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
}
}


int main(void){
double dt, t, tmax;
double X[4];
double Xderiv[4];

dt=0.01;
tmax=1000.0;

X[0]=1.0;
X[1]=0.0;
X[2]=0.0;
X[3]=-0.5;

ofstream fichier ("equadiffrk4cpp.out");

auto start = high_resolution_clock::now();//on commence a compter le temps de mesure

for(t=0.0;t<=tmax;t+=dt){
	rk4(t,X,dt,4,deriv);
	if((int)(round(t/dt))%10==0){//on n'ecrit qu'une valeur sur 10
	fichier <<t<<" "<<X[0]<<" "<<X[1]<<endl;
	}
}

auto stop = high_resolution_clock::now();//on arrete de compter le temps d'execution

fichier.close();

auto duration = duration_cast<microseconds>(stop - start); 
  
cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl; 
  
return 0; 

}

