#include<cmath>
#include<iostream>

using namespace std;

void funcx(double*, double*, double*, double&, double&, const double);
void funcy(double*, double*, double*, double&, double&, const double);

int main(){
  const double mu=0.012277471;
  double u[2]; //stores values of x and x_t
  double v[2]; //stores values of y and y_t
  double r, s;
  const double T=17.065216560157;
  double dt; 
  
  
  
  
  
  
  
  
  
  
  return0;
}

void funcx(double* f, double* u, double* v, double& r, double& s, const double mu){
  f[0]=u[1];
  f[1]=u[0]+2*v[1]-(1-mu)*(u[0]+mu)/(r*r*r)-mu*(u[0]-1+mu)/(s*s*s);
}
void funcy(double* f, double* u, double* v, double& r, double& s, const double mu){
  f[0]=v[1];
  f[1]=v[0]-2*u[1]-(1-mu)*v[0]/(r*r*r)-mu*v[0]/(s*s*s);
}