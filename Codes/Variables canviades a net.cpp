#include<iostream>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include <fstream>
#define num_cotxes 5
using namespace std;


int main(){// Es comença declarant alguns paràmetres
double C=20000.0;
double m=1500.0;
double L=4.0;
double d=26.0;
double norm_t=m*(L+d)/C;
// Es declaren els paràmetres en el cas on totes les variables seran diferents
//double c[num_cotxes]={21000,19000,20000,22000,18000};
//double M[num_cotxes]={1400.0,1950.0, 1165.0, 1280.0, 1100.0};
//double l[num_cotxes]={4.37,4.322,4.06, 4.227, 4.475};
//double D[num_cotxes]={20,26,26,24,18};
//int N[num_cotxes]={271,271,231,231,302};
int n=231;//es el mes gran de la llista N per tal de declarar les altres llistes amb aquest valor, en cas que N tinguis valors diferents
int N[num_cotxes]={n,n,n,n,n};// Es declaren les variables en el cas que tots els paràmetres siguin idèntics
double M[num_cotxes]={1500.0,1500.0, 1500.0, 1500.0, 1500.0};
double c[num_cotxes]={20000,20000,20000,20000,20000};
double D[num_cotxes]={26,26,26,26,26};
double l[num_cotxes]={4.0,4.0,4.0, 4.0, 4.0};
double L_norm=L/(L+d);
double l_norm[num_cotxes];
for(int i=0;i<num_cotxes;i++){
    l_norm[i]=l[i]/(L+d);// Es normalitza la distància
}
double d_min[num_cotxes-1];
double w=1.0;
double v_eq=120*norm_t*1000/((L+d)*3600);
double d_t=0.001;
double t_c=1/norm_t;
double t_final=16/norm_t;
int n_t=t_final/d_t;
double thau=n*d_t;
cout<<"thau="<<thau*norm_t<< "\n";
double K_1[num_cotxes-1];// Es defineixen les llistes que s'utilitzaran
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
double error[num_cotxes-1];
double error2=0;
double Ec_0[num_cotxes-1];
double W[num_cotxes-1];
for(int i=0;i<num_cotxes-1;i++){
	Ec_0[i]=M[i+1]*pow(v_eq,2)/(2*m);
	cout<< Ec_0[i] <<"\n";
	W[i]=0;
}
ofstream fitxer("Hola.txt");// Es creen 2 fitxers, un per imprimir les posicions i velocitats,  i l'altre l'error
ofstream fitxer2("Errors.txt");
for(int i=0;i<n;i++){
    x[i][0]=0;// s'imposen les condicions del cotxe que fa de condició de contorn
    v[i][0]=v_eq;
}
for(int i=0;i<=n;i++){
    for(int j=1;j<num_cotxes;j++){
    x[i][j]=x[i][j-1]-(l[j]+D[j])/(L+d);// s'imposa que els altres estiguin separats una certa distància
    v[i][j]=v_eq;
    }
}

fitxer << 0.000000 << "\t";
for(int j=0;j<num_cotxes;j++){
    fitxer<< x[0][j]*(L+d) << " \t " << v[0][j]*(L+d)/norm_t << "\t";// es comencen a imprimir variables
}
fitxer<< "\n";
for(double t=d_t;t<=(t_final+d_t);t+=d_t){
    t_c1+=d_t;
    fitxer<< t*norm_t<<"\t";
    if(t<=t_c){// s'inicialitza el procés abans de que es comenci a frenar
        for(int i=0;i<num_cotxes;i++){
         x[N[i]][i]=x[N[i]-1][i]+v_eq*d_t;
         v[N[i]][i]=v_eq;
         fitxer<< x[N[i]][i]*(L+d) << " \t " << v[N[i]][i]*(L+d)/norm_t << "\t";
    }
    for(int j=0;j<num_cotxes;j++){
            for(int i=1;i<N[j]+1;i++){
                v[i-1][j]=v[i][j];//el que es fa és dir que la posició del temps i-1 canvia a la delt temps i
                x[i-1][j]=x[i][j];
            }
        }
    fitxer<< "\n";
    }
    else{
    	// Es declaren tots els casos, tant la posició com la velocitat
        v[N[0]][0]=0.2*v_eq; //cas 1
        x[N[0]][0]=x[N[0]-1][0]+0.2*v_eq*d_t;//cas 1
        //v[N[0]][0]=v_eq*(1-(t-t_c)*norm_t*exp(1-(t-t_c)*norm_t));//cas 2
        //x[N[0]][0]=x_final+v_eq*((t-t_c+1/norm_t)*exp((t_c-t+1/norm_t)*norm_t)-t_c+t-exp(1)/norm_t);//cas 2
        //v[N[0]][0]=v_eq*(1-0.8*pow(sin(w*(t-t_c)),2)); //cas 3
        //x[N[0]][0]=x_final+v_eq*((t-t_c)-0.2*(-sin((-t_c+t)*2*w)+2*w*(t-t_c))/w);//cas 3
        fitxer<< x[N[0]][0]*(L+d) << "\t " << v[N[0]][0]*(L+d)/norm_t << "\t";
        for(int j=0;j<num_cotxes-1;j++){
        		// s'aplica el mètode de Runge-Kutta
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
                fitxer<< x[N[j+1]][j+1]*(L+d) << " \t " << v[N[j+1]][j+1]*(L+d)/norm_t << "\t";
        }
        for(int j=0;j<num_cotxes;j++){
            for(int i=1;i<N[j]+1;i++){
                v[i-1][j]=v[i][j];
                x[i-1][j]=x[i][j];
            }
    	}
        for(int i=0;i<num_cotxes-1;i++){// Es calcula l'error i s'imprimeix
        	W[i]+=M[i+1]*K_1[i]*d_t*(L_1[i]+2*L_2[i]+2*L_3[i]+L_4[i])/(6*m);
        	error[i]=100*(Ec_0[i]+W[i]-M[i+1]*pow(v[N[i+1]][i+1],2)/(2*m))/(Ec_0[i]);
        	fitxer2 << error[i] << "\t";
		}
		for(int i=1;i<num_cotxes;i++){// es crea un if que ens surti del bucle temporal si els cotxes xoquen
            if((d+L)*(x[N[i-1]][i-1]-x[N[i]][i])<l[i-1]){
                cout << "Els cotxes"<< "\t" << i-1 << "\t " << "i" << " \t"<< i <<"\t"<<" s'han xocat i tot ha parat en el temps:"<< t_c1*norm_t;
				goto fora;
			}
		}

    fitxer<< "\n";
    fitxer2 << t*norm_t << "\n";
    }
}
fora:
fitxer.close();
fitxer2.close();
return 0;
}

