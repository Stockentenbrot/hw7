#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

void f(double t, double *y,double *k)
{
  const double mu=0.012277471;
  double einsdurchrhochdrei = 1.0/pow(sqrt((y[0]+mu)*(y[0]+mu)+y[1]*y[1]),3);
  double einsdurchshochdrei = 1.0/pow(sqrt((y[0]-1.0+mu)*(y[0]-1.0+mu)+y[1]*y[1]),3);
  k[0] = y[2];
  k[1] = y[3];
  k[2] = y[0]+2.0*y[3]-(1.0-mu)*(y[0]+mu)*einsdurchrhochdrei-mu*(y[0]-1.0+mu)*einsdurchshochdrei;
  k[3] = y[1]-2.0*y[2]-(1.0-mu)*y[1]*einsdurchrhochdrei-mu*y[1]*einsdurchshochdrei;
  return;
}

void Berechne_k(const int dim,const double t,const double dt,const int s, double **k, double **a, double *c, double *y) /// Indizes wie in der Formel bei Wikipedia ;)
{
  int j,l,d;
  double Summe;
  double *ytemp;
  
  ytemp = (double*) malloc(dim*sizeof(double));
  for (j=0; j < s; j++) /// alle k
  {
    for (d=0; d < dim; d++) /// in allen Dimensionen
    {
      Summe = 0.0;
      for (l=0; l < s; l++)
      {
	Summe += a[j][l]*k[l][d];
      }
      ytemp[d] = y[d] + dt*Summe;
      
    }
    f(t + c[j]*dt,ytemp,k[j]);
  }
  free(ytemp);
  return;
}

void Berechne_y(const int dim,const double dt,const int s, double **k, double *b, double *y,double *yneu)
{
  double Summe;
  int d,j;
  
  for (d=0; d < dim; d++) /// in allen Dimensionen
  {
      Summe = 0.0;
      for (j=0; j < s; j++)
      {
	Summe += b[j]*k[j][d];
      }
      yneu[d] = y[d] + dt*Summe;
  }
  return;
}

void Step(const int dim,const double t,const double dt,const int s, double **a,double *b, double *c, double *y, double *yneu, double **k)
{
  Berechne_k(dim,t,dt,s,k,a,c,y);
  Berechne_y(dim,dt,s,k,b,y,yneu);
  return;
}

int main(void)
{
double dt=1E-3,Tol=1E-5;
const int s = 7;
const int dim = 4;
double t=0.0;
const double L=17.065216560157;

double temp;
double max = 0.0;
int i,j;


/// Plain C-style ;) denn "a = new double[s][s];" wollte einfach nicht funktionieren...
double **a;
a = (double**)malloc(s*sizeof(double*));
for (i = 0; i < s; i++)
{
  a[i] = (double*)malloc(s*sizeof(double));
}

double **k;
k = (double**)malloc(s*sizeof(double*));
for (i = 0; i < s; i++)
{
  k[i] = (double*)malloc(dim*sizeof(double));
}

for (i=0; i < s; i++)
{
  for (j=0; j < s; j++)
  {
    a[i][j] = 0.0;
  }
}
a[1][0] = 1.0/5.0;
a[2][0] = 3.0/40.0; a[2][1] = 9.0/40.0;
a[3][0] = 44.0/45.0;  a[3][1] = -56.0/15.0;  a[3][2] = 32.0/9.0;
a[4][0] = 19372.0/6561.0;  a[4][1] = -25360.0/2187.0;  a[4][2] = 64448.0/6561.0; a[4][3] = -212.0/729.0;
a[5][0] = 9017.0/3168.0;  a[5][1] = -355.0/33.0;  a[5][2] = 46732.0/5247.0; a[5][3] = 49.0/176.0; a[5][4] = -5103.0/18656.0;
a[6][0] = 35.0/384.0;  a[6][1] = 0.0;  a[6][2] = 500.0/1113.0; a[6][3] = 125.0/192.0; a[6][4] = -2187.0/6784.0; a[6][5] = 11.0/84.0;

double b4[s];
b4[0] = 35.0/384.0; b4[1] = 0.0; b4[2] = 500.0/1113.0; b4[3] = 125.0/192.0; b4[4] = -2187.0/6784.0; b4[5] = 11.0/84.0; b4[6] = 0.0;
double b5[s];
b5[0] = 5179.0/57600.0; b5[1] = 0.0; b5[2] = 7571.0/16695.0; b5[3] = 393.0/640.0; b5[4] = -92097.0/339200.0; b5[5] = 187.0/2100.0; b5[6] = 1.0/40.0;
double c[s];
c[0] = 0.0; c[1] = 1.0/5.0; c[2] = 3.0/10.0; c[3] = 4.0/5.0; c[4] = 8.0/9.0; c[5] = 1.0; c[6] = 1.0;

double y[dim];
double yneu4[dim];
double yneu5[dim];
y[0]=0.994; y[1]=0.0; y[2]=0.0; y[3]=-2.00158510637908;

ofstream Datei("Ergebnis.dat");
Datei << t << '\t' << dt <<'\t' << y[0] << '\t' << y[1] << '\t' << endl;

while(t < L)
{
  Step(dim,t,dt,6,a,b4,c,y,yneu4,k);
  Step(dim,t,dt,7,a,b5,c,y,yneu5,k);
  t +=dt;
  for(i=0; i<dim; i++)
  {
    y[i] = yneu4[i];
  }
  Datei << t << '\t' << dt <<'\t' << y[0] << '\t' << y[1] << '\t' << endl;

  max = 0.0;
  for(i=0; i<dim; i++)
  {
    temp = abs(yneu4[i]-yneu5[i]);
    if (temp>max)
    {
      max = temp;
    }
  }
  dt *= pow(Tol/max, 0.2); 
}
Datei.close();


for (i = 0; i < s; i++)
{
  free(k[i]);
}
free(k);


for (i = 0; i < s; i++)
{
  free(a[i]);
}
free(a);

exit(0);
}




