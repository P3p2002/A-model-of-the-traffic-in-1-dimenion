#include<iostream>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include <fstream>
#define num_cotxes 5
using namespace std;

int main(){
// Es decalren diverses variables
double C=20000.0;
double m=1500.0;
double L=4.0;
double d=26.0;
int n=150;//es el mes gran de la llista N per tal de declarar les altres llistes amb aquest valor, si és que N té valors diferents
// Es creen aquestes llistes per utilitzar amb les variables diferents
int N[num_cotxes]={50,100,150,150,150};
double M[num_cotxes]={1400.0,1950.0, 1165.0, 1280.0, 1100.0};
double c[num_cotxes]={21000,19000,20000,22000,18000};
double l[num_cotxes]={4.37,4.322,4.06, 4.227, 4.475};
double D[num_cotxes]={20,26,26,24,18};
// Es creen aquestes llistes per utilitzar quan les variables siguin iguales
//double D[num_cotxes]={26,26,26,26,26};
//double l[num_cotxes]={4.0,4.0,4.0, 4.0, 4.0};
//double c[num_cotxes]={20000,20000,20000,20000,20000};
//double M[num_cotxes]={1500.0,1500.0, 1500.0, 1500.0, 1500.0};
//int N[num_cotxes]={n,n,n,n,n};
// Es declaren diverses variables
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
double t_final=16/norm_t;
int n_t=t_final/d_t;
double thau=n*d_t;
double t_parada=1000/norm_t;//ho estableixo amb un valor exageradament gran per tal de que mes endavanat en el codi m'ho pari quan ho vull
double t_arrancada=1.5/norm_t;
double x_parada;
printf("thau=%lf\n",thau*norm_t);
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
double t_c2;
ofstream fitxer("Hola.txt");
// Es posen els valors del cotxe que fa de condició de contorn
for(int i=0;i<n;i++){
    x[i][0]=0;
    v[i][0]=v_eq;
}
// Imposo que els cotxes estiguin a una certa distància
for(int i=0;i<n;i++){
    for(int j=1;j<num_cotxes;j++){
    x[i][j]=x[i][j-1]-(l[j]+D[j])/(L+d);
    v[i][j]=v_eq;
    }
}
fitxer<< 0.000000 <<"\t";
for(int j=0;j<num_cotxes;j++){
    fitxer<< x[0][j]*(L+d) << " \t " << v[0][j]*(L+d)/norm_t << "\t";
}
fitxer <<"\n";
for(double t=d_t;t<=(t_final+d_t);t+=d_t){
    t_c1+=d_t;
    fitxer << t*norm_t << "\t";
    if(t<=t_c){// Començo inicialitzant tot abans de que comenci a frenar el primer cotxe
        for(int i=0;i<num_cotxes;i++){
         x[N[i]][i]=x[N[i]-1][i]+v_eq*d_t;
         v[N[i]][i]=v_eq;
         fitxer << x[N[i]][i]*(L+d) << "\t" << v[N[i]][i]*(L+d)/norm_t << "\t";
    }
    for(int j=0;j<num_cotxes;j++){
            for(int i=1;i<N[j]+1;i++){
                v[i-1][j]=v[i][j];//el que es fa és dir que la posició del temps i-1 canvia a la delt temps i
                x[i-1][j]=x[i][j];
            }
        }
    fitxer <<"\n";
    }
    else{
        if(v[N[0]][0]>=0.01 && t<t_parada){//el fet d'introduir aqui un t_parada molt gran em fara que mentres la velocitat no sigui nula em faci aquest bucle, pero un cop es faci 0 no em torni en aquest bucle
            // Aquí estic imposant que el cotxi freni i es quedi proper a una velocitat nul·la
			v[N[0]][0]=v_eq*(1-(t-t_c)*norm_t*exp(1-(t-t_c)*norm_t));
            x[N[0]][0]=x_final+v_eq*((t-t_c+1/norm_t)*exp((t_c-t+1/norm_t)*norm_t)-t_c+t-exp(1)/norm_t);
            if(v[N[0]][0]<0.01){
                t_parada=t;// em guardo el temps de parada per tal de que no es torni a repetir el bucle
                t_c2=t_parada+t_arrancada;
                x_parada=x[N[0]][0];// em guardo la posició on ha parat el cotxe que fa de condició de contorn per més endavant
            }
        }
        else if(t<t_c2){// Imposo que estigui frenat el temps que vull el primer cotxe
            v[N[0]][0]=0;
            x[N[0]][0]=x_parada;
        }
        else if(t_c2<=t){// Després de que estigui frenat imposo que arrenqui tal com li he demanat
            v[N[0]][0]=v_eq*log(1+(t-t_c2)*norm_t)/4;
            x[N[0]][0]=x_parada+v_eq*(log(1+(t-t_c2)*norm_t)*(1+(t-t_c2)*norm_t)+(t_c2-t)*norm_t)/(4*norm_t);
        }
        else if(fabs(v[N[0]][0]-v_eq)<0.001 && t_c2<=t){ // Impooso que la velocitat es mantingui constant un cop ha arribat a una prou alta
            v[N[0]][0]=v_eq;
            x[N[0]][0]=x[N[0]-1][0]+v_eq*d_t;
        }
        fitxer << x[N[0]][0]*(L+d) << "\t" << v[N[0]][0]*(L+d)/norm_t << "\t";
        for(int j=0;j<num_cotxes-1;j++){
        	    // Aplico Runge-Kutta 
                K_1[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]-x[0][j])-l_norm[j]));
                L_1[j]=v[N[j+1]-1][j+1];
                K_2[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_1[j]*d_t/2-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]-x[0][j])+L_1[j]*d_t/2-l_norm[j]));
                L_2[j]=v[N[j+1]-1][j+1]+K_1[j]*d_t/2;
                K_3[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_2[j]*d_t/2-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]-x[0][j])+L_2[j]*d_t/2-l_norm[j]));
                L_3[j]=v[N[j+1]-1][j+1]+K_2[j]*d_t/2;
                K_4[j]=-1*(m*c[j+1])/(M[j+1]*C)*(v[N[j+1]-1][j+1]+K_3[j]*d_t-v[0][j])/(fabs(fabs(x[N[j+1]-1][j+1]-x[0][j])+L_3[j]*d_t-l_norm[j]));
                L_4[j]=v[N[j+1]-1][j+1]+K_3[j]*d_t;
                v[N[j+1]][j+1]=v[N[j+1]-1][j+1]+d_t*(K_1[j]+2*K_2[j]+2*K_3[j]+K_4[j])/6;
                x[N[j+1]][j+1]=x[N[j+1]-1][j+1]+d_t*(L_1[j]+2*L_2[j]+2*L_3[j]+L_4[j])/6;
                fitxer<< x[N[j+1]][j+1]*(L+d) << "\t" << v[N[j+1]][j+1]*(L+d)/norm_t << "\t";
        }
        for(int j=0;j<num_cotxes;j++){
            for(int i=1;i<N[j]+1;i++){
                v[i-1][j]=v[i][j];
                x[i-1][j]=x[i][j];
            }
        }
        for(int j=0;j<num_cotxes;j++){// Aquí creo un bucle per parar la resta de cotxes
            if((L+d)*v[N[j]][j]/norm_t<0.01){
                v[N[j]][j]=0;
            }
        }
        for(int i=1;i<num_cotxes;i++){// Creo un if que em pari el que estic fent en cas de que xoquin els cotxes
            if((d+L)*(x[N[i-1]][i-1]-x[N[i]][i])<l[i-1]){
                cout<<"Els cotxes"<< i << "i" << i <<" s'han xocat i tot ha parat en el temps:" << t_c1*norm_t;
                goto fora;
            }

        }
    fitxer<<"\n";
    }
}
fora:
fitxer.close();
return 0;
}

