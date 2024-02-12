#include<iostream>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include <fstream>
#define num_cotxes 5
using namespace std;

int main(){
double C=20000.0;// es comencen declarant variables
double m=1500.0;
double L=4.0;
double d=26.0;
int n;
// aquestes llistes es fan servir en cas que els paràmetres siguin idèntics
double M[num_cotxes]={1500.0,1500.0, 1500.0, 1500.0, 1500.0};
double c[num_cotxes]={20000,20000,20000,20000,20000};
double l[num_cotxes]={4.0,4.0,4.0, 4.0, 4.0};
double D[num_cotxes]={26,26,26,26,26};
// aquestes llistes es fan servir en cas que els paràmetres siguin diferents
//double D[num_cotxes]={20,26,26,24,18};
//double l[num_cotxes]={4.37,4.322,4.06, 4.227, 4.475};
//double c[num_cotxes]={21000,19000,20000,22000,18000};
//double M[num_cotxes]={1400.0,1950.0, 1165.0, 1280.0, 1100.0};
// es defineixen diverses variables
double L_norm=L/(L+d);
double l_norm[num_cotxes];
for(int i=0;i<num_cotxes;i++){
    l_norm[i]=l[i]/(L+d);
}
double norm_t=m*(L+d)/C;
double w=1.0;
double v_eq=120*norm_t*1000/((L+d)*3600);
double d_t=0.001;
double t_c=1/norm_t;
double t_final=22/norm_t;
int n_t=t_final/d_t;
double tau_max=0.45;//Esta normalitzat per tal de que sigui un multiple de d_t
int n_0=5;
double T_inicial=n_0*d_t;
int n_tau=tau_max/d_t-n_0;
double d_min[n_tau][num_cotxes-1];
double TAU[n_tau];
double min_dmin[n_tau];
TAU[0]=T_inicial;
int t_tau=0;
// Es crea un valor de la distància mínim que ja sabem que serà major a la distància mínima real
for(int i=0;i<=n_tau;i++){
    for(int j=0;j<num_cotxes-1;j++){
        d_min[i][j]=(l[j]+D[j])/(L+d);
    }
}
for(double T=T_inicial;T<tau_max;T+=d_t){// es crea el bucle dels temps de reacció
    n=T/d_t;// és el n que farà de multiple del temps de reacció
    // es creen  les llistes que ens interesen cada cop que es fa el bucle del temps de reacció
    int N[num_cotxes]={n,n,n,n,n};
    double K_1[num_cotxes-1];
    double L_1[num_cotxes-1];
    double K_2[num_cotxes-1];
    double L_2[num_cotxes-1];
    double K_3[num_cotxes-1];
    double L_3[num_cotxes-1];
    double K_4[num_cotxes-1];
    double L_4[num_cotxes-1];
    double x[n+1][num_cotxes];
    double x_final=v_eq/norm_t;
    double v[n+1][num_cotxes];
    double t_c1=0;
    for(int i=0;i<n;i++){// s'inicialitza el cotxe que fa de condició de contorn
        x[i][0]=0;
        v[i][0]=v_eq;
    }
    for(int i=0;i<n;i++){
        for(int j=1;j<num_cotxes;j++){
        x[i][j]=x[i][j-1]-(l[j]+D[j])/(L+d);// es fa que els cotxes estiguin separats una certa distància
        v[i][j]=v_eq;
        }
    }
    for(int j=0;j<num_cotxes;j++){
    }
    for(double t=d_t;t<=(t_final+d_t);t+=d_t){// es crea el bucle temporal ara
        t_c1+=d_t;
        if(t<=t_c){
            for(int i=0;i<num_cotxes;i++){// es creen els valors abans de que comenci a frenar el primer cotxe
             x[N[i]][i]=x[N[i]-1][i]+v_eq*d_t;
             v[N[i]][i]=v_eq;
        }
        for(int j=0;j<num_cotxes;j++){
                for(int i=1;i<N[j]+1;i++){
                    v[i-1][j]=v[i][j];//el que es fa és dir que la posició del temps i-1 canvia a la delt temps i
                    x[i-1][j]=x[i][j];
                }
            }
        }
        else{
        	// Es declaren tots els cassos, tant de la posició com de la velocitat
            v[N[0]][0]=0.2*v_eq; //cas 1
            x[N[0]][0]=x[N[0]-1][0]+0.2*v_eq*d_t;//cas 1
            //v[N[0]][0]=v_eq*(1-(t-t_c)*norm_t*exp(1-(t-t_c)*norm_t));//cas 2
            //x[N[0]][0]=x_final+v_eq*((t-t_c+1/norm_t)*exp((t_c-t+1/norm_t)*norm_t)-t_c+t-exp(1)/norm_t);//cas 2
            //v[N[0]][0]=v_eq*(1-0.8*pow(sin(w*(t-t_c)),2)); //cas 3
        	//x[N[0]][0]=x_final+v_eq*((t-t_c)-0.2*(-sin((-t_c+t)*2*w)+2*w*(t-t_c))/w);//cas 3
            for(int j=0;j<num_cotxes-1;j++){
            	// S'aplica Runge-Kutta
                    K_1[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]-x[0][j])-l_norm[j]));
                    L_1[j]=v[N[j+1]-1][j+1];
                    K_2[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_1[j]*d_t/2-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]+L_1[j]*d_t/2-x[0][j])-l_norm[j]));
                    L_2[j]=v[N[j+1]-1][j+1]+K_1[j]*d_t/2;
                    K_3[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_2[j]*d_t/2-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]+L_2[j]*d_t/2-x[0][j])-l_norm[j]));
                    L_3[j]=v[N[j+1]-1][j+1]+K_2[j]*d_t/2;
                    K_4[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_3[j]*d_t-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]+L_3[j]*d_t-x[0][j])-l_norm[j]));
                    L_4[j]=v[N[j+1]-1][j+1]+K_3[j]*d_t;
                    v[N[j+1]][j+1]=v[N[j+1]-1][j+1]+d_t*(K_1[j]+2*K_2[j]+2*K_3[j]+K_4[j])/6;
                    x[N[j+1]][j+1]=x[N[j+1]-1][j+1]+d_t*(L_1[j]+2*L_2[j]+2*L_3[j]+L_4[j])/6;
            }
            for(int j=0;j<num_cotxes;j++){
                for(int i=1;i<N[j]+1;i++){
                    v[i-1][j]=v[i][j];
                    x[i-1][j]=x[i][j];
                }
            }
            for(int i=0;i<num_cotxes-1;i++){// Imposo la condició que si les distàncies són menor a la mínima que ja teniem la guarda
                if(fabs(x[N[i+1]][i+1]-x[N[i]][i]+l_norm[i])<d_min[t_tau][i]){
                    d_min[t_tau][i]=fabs(x[N[i+1]][i+1]-x[N[i]][i]+l_norm[i]);
                }
            }
            for(int i=1;i<num_cotxes;i++){// Si els cotxes xoquen surto, però torno al bucle del temps de reacció
            	if((d+L)*(x[N[i-1]][i-1]-x[N[i]][i])<l[i-1]){
                	cout << "Els cotxes"<< "\t" << i-1 << "\t " << "i" << " \t"<< i <<"\t"<<" s'han xocat i tot ha parat en el temps:"<< t_c1*norm_t;
                	t_tau+=1;
					goto fora;
            	}

        	}
        }
    }
    t_tau+=1;
    TAU[t_tau]=T;
}
fora:
min_dmin[t_tau]=0;
// Agafo el mínim de les distàncies mínimes
ofstream fitxer("Distancies minimes.txt");
for(int i=1;i<t_tau;i++){
    min_dmin[i]=d_min[i][0];
    for(int j=0;j<num_cotxes-1;j++){
        if(d_min[i][j]<min_dmin[i]){
            min_dmin[i]=d_min[i][j];
        }
    }
    fitxer<< TAU[i]*norm_t << " \t " << min_dmin[i]*(L+d)<< "\n";
}
fitxer.close();
return 0;
}
