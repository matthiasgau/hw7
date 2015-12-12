#include<cmath>
#include<iostream>
#include<fstream>

using namespace std;

void func(double*, double, double, double, double, const double);
void k_fill(double*, double*, double*, double*, double*, double*, double*, double*, double, const double);
void Kutta4(double*, double*, double*, double*, double*, double*, double*, double*, double, const double);
void Kutta5(double*, double*, double*, double*, double*, double*, double*, double*, double, const double);
void step_size(double*, double*, double*, double&, double&, double&);

int main(){
  const double mu=0.012277471;
  double u[4]; //stores values of x, x_t, y and y_t for RK4
  double v[4]; //stores the same for Rk5
  const double T=17.065216560157;
  double dt; 
  double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4];
  double norm[4]; //stores the norm of x, x_t, y and y_t with respect to both RK methods 
  double max1, max2, maximum; //we need to store to compare the norm values
  double t=0.0; //starting time
  u[0]=0.994;u[1]=0;u[2]=0;u[3]=-2.00158510637908;
  v[0]=0.994;v[1]=0;v[2]=0;v[3]=-2.00158510637908;
  double w[4]; //need to store values of RK4 for the embedded RK
  double x[4]; //storage of Rk5 results
  
  double Tol=1e-5; //toleranz may be changed
  double q=0.7; //safty factor
  
  
  ofstream peter ("Kutter.txt");
  
  peter<< t << "\t" << 0 << "\t" << u[0] << "\t" << u[2] <<endl;
  
  while(t<T){
    dt=Tol; //start value of dt
    for(int i=0;i<4;i++) //store the old values of u
      w[i]=u[i];
    for(int i=0;i<4;i++) // and v
      x[i]=v[i];
    Kutta4(k1,k2,k3,k4,k5,k6,k7,u,dt,mu); //now u hosts the new values given by RK4
    Kutta5(k1,k2,k3,k4,k5,k6,k7,v,dt,mu); //v hosts values given by RK5
    step_size(norm, u,v,max1,max2,maximum); //calculates the maximimum difference of both RK methods for the old dt
    dt*=pow((Tol/maximum),0.2); //calculates new step size
    for(int i=0;i<4;i++) //need the old values of u and v 
      u[i]=w[i];         //since we want to calculate just one step with the 
    for(int i=0;i<4;i++) //fitted step size
      v[i]=x[i];
    Kutta4(k1,k2,k3,k4,k5,k6,k7,u,dt,mu);
    Kutta5(k1,k2,k3,k4,k5,k6,k7,v,dt,mu);
    t+=dt;
    
    peter<< t << "\t" << dt << "\t" << u[0] << "\t" << u[2] <<endl;
  }
   
  peter.close();
  
  
  return 0;
}

void func(double* f, double u0, double u1, double u2, double u3, const double mu){
  f[0]=u1;
  f[1]=u0+2.0*u3-(1.0-mu)*(u0+mu)/pow(sqrt(pow(u0+mu,2)+pow(u2,2)),3)-mu*(u0-1.0+mu)/(pow(sqrt(pow(u0-1.0+mu,2)+pow(u2,2)),3));
  f[2]=u3;
  f[3]=u2-2.0*u1-(1.0-mu)*u2/(pow(sqrt(pow(u0+mu,2)+pow(u2,2)),3))-mu*u2/(pow(sqrt(pow(u0-1.0+mu,2)+pow(u2,2)),3));
}
void k_fill(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* u, double dt, const double mu){
  double u0=u[0], u1=u[1], u2=u[2], u3=u[3];
  const double a21=1.0/5.0;
  const double a31=3.0/40.0, a32=9.0/40.0;
  const double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
  const double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
  const double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0;
  const double a71=35.0/385.0, a72=0.0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;
  func(k1, u0, u1, u2, u3, mu);
  func(k2, u0+dt*a21*k1[0], u1+dt*a21*k1[1], u2+dt*a21*k1[2], u3+dt*a21*k1[3], mu);
  func(k3, u0+dt*(a31*k1[0]+a32*k2[0]), u1+dt*(a31*k1[1]+a32*k2[1]), u2+dt*(a31*k1[2]+a32*k2[2]), u3+dt*(a31*k1[3]+a32*k2[3]), mu);
  func(k4, u0+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]), u1+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]), u2+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), u3+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]), mu);
  func(k5, u0+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]), u1+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]), u2+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), u3+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]), mu);
  func(k6, u0+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]), u1+dt*(a61*k1[1]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]), u2+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), u3+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]), mu);
  func(k7, u0+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]), u1+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]), u2+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), u3+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]), mu);
}
void Kutta4(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* u, double dt, const double mu){  
  const double b1=5179.0/57600.0, b2=0.0, b3=7571.0/16695.0, b4=393.0/640.0, b5=-92097.0/339200.0, b6=187.0/2100.0, b7=1.0/40.0;
  k_fill(k1, k2, k3, k4, k5, k6, k7, u, dt, mu);
  for(int i=0;i<4;i++)
    u[i]+=dt*(b1*k1[i]+b2*k2[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]+b7*k7[i]);
}
void Kutta5(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* u, double dt, const double mu){
  const double b1=35.0/384.0, b2=0.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0, b7=0.0;
  k_fill(k1, k2, k3, k4, k5, k6, k7, u, dt, mu);
  for(int i=0;i<4;i++)
    u[i]+=dt*(b1*k1[i]+b2*k2[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]+b7*k7[i]);
}
void step_size(double* norm, double* u, double* v, double& max1, double& max2, double& maximum){
  for(int i=0;i<4;i++)
    norm[i]=abs(u[i]-v[i]);
  max1=max(norm[0],norm[1]);
  max2=max(norm[2],norm[3]);
  maximum=max(max1,max2);
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  