/* Dieses Script rechnet die Punkte für die Relation M(R) für gegebenes P_c aus.*/



#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

// dp/dr
double func1(double r, double p, double m); // hier werden die Funktionen deklariert
// dm/dr
double func2(double r, double p); 


void funcmr(double p0);

double const c_g = 6.6743*pow(10,-11);  // Gravitationskonstante G
//double pn = 1.164*pow(10,31);           // Dies ist der Anfangswert des Drucks P(0)
double mn = 0.0;                       // Der Anfangswert der eingeschlossenen Masse m(0)
double rn = pow(10,-30);                       // Der Startwert des Radius
double const c_delr = 0.01;                // Die Schrittweite delta r
double const c_c = 299792458.0;         // Lichtgeschwindigkeit c
double const c_gamma = 1.0/(1+1/1.0);   // Der Exponent des Polytropen EOS mit dem Polytropenindex 0.5  rho = (P/K)^gamma    
double const c_k = 5*pow(10,-3);                 // Konstante K
double const c_pi= 3.14159265359;                // pi






// hier werden die Koeffizienten des Runge-Kutta-Verfahrens deklariert
double K1f; // f steht hier für die DGL dp/dr
double K1g; // g steht hier für die DGL dm/dr 
double K2f;
double K2g;
double K3f;
double K3g;
double K4f;
double K4g;





    // hier wird die Funktion dp/dr definiert, welche für (r,p,m) einen Wert ausgibt
    double func1(double r, double p, double m){
    double F1 = - (c_g*(m + 4* c_pi*pow(r,3)*p/(pow(c_c,2)))*((pow((p/c_k),c_gamma))+p/pow(c_c,2)))/(r*(r-(2*c_g*m)/pow(c_c,2)));
    return F1;
    }

    double func2(double r, double p){
    double F2 = 4*c_pi*pow(p/c_k,c_gamma)*pow(r,2);
    return F2;
    }
double M, R;

//Diese Funktion weist M und R ihre Werte für ein gegebens P_c zu.
void funcmr(double pn){

    rn=pow(10,-30);
    mn=0.0;
    while (pn > 1){
        
        K1f = c_delr*func1(rn, pn, mn); //Die K_i werden hier berechnet
        K1g = c_delr*func2(rn,pn);
        K2f = c_delr*func1(rn+ (1.0/2) *c_delr,pn+(1.0/2)*K1f,mn+(1.0/2)*K1g);
        K2g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K1f);
        K3f = c_delr*func1(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f,mn+(1.0/2)*K2g);
        K3g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f);
        K4f = c_delr*func1(rn+c_delr,pn+K3f,mn+K3g);
        K4g = c_delr*func2(rn+c_delr,pn+K3f);
        
        M= mn;
        R= rn;
        pn += (1.0/6) *K1f+(1.0/3)*K2f+(1.0/3)*K3f+(1.0/6)*K4f; //Der Wert p(r+delr) ergibt sich hier
        mn += (1.0/6) *K1g+(1.0/3)*K2g+(1.0/3)*K3g+(1.0/6)*K4g; //Der Wert m(r+delr) ergibt sich hier
        rn +=c_delr;



     
        
      
    
    }
}


// In der Main fct wird funcmr() aufgerufen und die Datenpunkte werde abgespeichert.
int main(){
     ofstream file1;
    file1.open ("datenpunkte2.csv");
     for (double p = 1.164*pow(10,31) ; p < 7.562*pow(10,45); p = p*1.1)
     {
        funcmr(p);
        file1 <<R<<","<<M<<endl;
     }
     file1.close();

    return 0;
}
