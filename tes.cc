





#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <iomanip>
#include <cmath>
#include <array>
#include <ostream>
#include <fstream>

using namespace std;

namespace {
}

const int N=4;
const int Ncell=2;



typedef struct
{
	int type;
	double x[N];
	double x_prev[N];
	double size_prev;
        double size;
} Cell;


double dt=0.05;
int Tmax=1000000;
double c=0.2;
double D0=0.1;
double Venv=1.0;
double Dcell=0.1;
double k01=0.1;
double k12=0.1;
double k23=0.1;
const double K[Ncell][N][N]={{{-k01,0,0,0},{k01,-k12,0,0},{0,k12,-k23,0},{0,0,k23,0}},{{-10.0*k01,0,0,0},{10.0*k01,-10.0*k12,0,0},{0,10.0*k12,-10.0*k23,0},{0,0,10.0*k23,0}}};
const int Cat[Ncell][N][N]={{{2,0,0,0},{2,3,0,0},{0,3,1,0},{0,0,1,0}},{{2,0,0,0},{2,3,0,0},{0,3,1,0},{0,0,1,0}}};
const double sigma[Ncell][N]={{Dcell,0,0,0},{Dcell,0,0,0}};


void init(Cell *cell,double *x0_prev){
  for(int i=0;i<Ncell;++i){
    cell[i].size_prev=1.0/(1.0*Ncell);
    for(int j=0;j<N;++j){cell[i].x[j]=1.0/(1.0*N);cell[i].x_prev[j]=1.0/(1.0*N);}

  }
  for(int j=0;j<N;++j){x0_prev[j]=0;}
  x0_prev[0]=c;
}

void react_cell(Cell *cell, double *x0_prev){
  for(int i=0;i<Ncell;++i){
    double lambda=0;
    for(int j =0;j<N;++j){lambda+=(sigma[i][j]*(x0_prev[j]-cell[i].x_prev[j]));}
    for(int j =0;j<N;++j){
      cell[i].x[j]=cell[i].x_prev[j]+dt*(sigma[i][j]*(x0_prev[j]-cell[i].x_prev[j]))-dt*cell[i].x_prev[j]*lambda;
      for(int k =0;k<N;++k){cell[i].x[j]+=dt*K[i][j][k]*cell[i].x[k]*cell[i].x[Cat[i][j][k]];}
    }
  }
}

void react_env(Cell *cell,double *x0, double *x0_prev){
  for(int i=0;i<N;++i){x0[i]=x0_prev[i];}
  x0[0]=x0_prev[0]+D0*(c-x0_prev[0])*dt;
  for(int i=0;i<N;++i){
    for(int j=0;j<Ncell;++j){
      x0[i]+=sigma[j][i]*(cell[j].x_prev[i]-x0_prev[i])*cell[j].size_prev*dt/Venv;
    }
  }
}

void growth(Cell *cell, double *x0_prev){
  double lambda[Ncell];
  for(int i =0;i<Ncell;++i){
    lambda[i]=0;
    for(int j =0;j<N;++j){lambda[i]+=(sigma[i][j]*(x0_prev[j]-cell[i].x_prev[j]));}
  }
  double Gamma=0;
  for(int i =0;i<Ncell;++i){Gamma+=lambda[i]*cell[i].size_prev;}
  for(int i =0;i<Ncell;++i){
    cell[i].size=cell[i].size_prev+dt*(cell[i].size_prev*lambda[i]-cell[i].size_prev*Gamma);}
  }
void update(Cell *cell,double *x0, double *x0_prev){
  for(int i=0;i<Ncell;++i){
    for(int j=0;j<N;++j){cell[i].x_prev[j]=cell[i].x[j];}
    cell[i].size_prev=cell[i].size;
  }
  for(int j=0;j<N;++j){x0_prev[j]=x0[j];}
}

  int main(){
    double x0[N];
    double x0_prev[N];

    Cell cell[Ncell];
    init(cell,x0_prev);
    for(int t=0;t<Tmax;++t){
      react_cell(cell,x0_prev);
      react_env(cell,x0,x0_prev);
      growth(cell,x0_prev);
      update(cell,x0,x0_prev);



      if(t%1000==0){
      cout<<t*dt<<" ";
      for(int i=0;i<Ncell;++i){
	for(int j=0;j<N;++j){cout<<cell[i].x[j]<<" ";}
      }
      for(int j=0;j<N;++j){cout<<x0[j]<<" ";}
      for(int i=0;i<Ncell;++i){cout<<cell[i].size<<" ";}
      cout<<endl;
      // gnuplot出力の際のindex
      //1: 時間
      //2+i*N+j: i番目のcellのj番目のchemicalの濃度。　
      //2+N*Ncell+j: j番目のchemicalの外部濃度。　
      //2+N*(Ncell+1)+i:　i番目のcellのsize
      // ex) plot "..." u 1:14 w l, "..." u 1:15 w l なら細胞1と2の成長率をplotする。...は出力ファイル名
    }
    }
}




