/* In diesem script werden die TOV-Gleichungen numerisch gelöst (Das Skript diente nur dazu zu schauen, ob 
   sich etwas sinnvolles ergibt, weshalbt es noch relativ unsauber ist).
Akutueller Stand:   Der Druck wird kleiner, fällt jedoch zu langsam, sodass der Radius zu sehr steigt. 
                    Der Druck fällt schlagartig und die Konsole gibt für den Druck nan als letzten Wert aus. 
 */



#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

// dp/dr
double func1(double r, double p, double m); // hier werden die Funktionen deklariert
// dm/dr
double func2(double r, double p); 

double const c_g = 6.6743*pow(10,-11);  // Gravitationskonstante G
double pn = 1.164*pow(10,31);           // Dies ist der Anfangswert des Drucks P(0)
double mn = 0.0;                       // Der Anfangswert der eingeschlossenen Masse m(0)
double rn = pow(10,-30);                       // Der Startwert des Radius
double const c_delr = 1;                // Die Schrittweite delta r
double const c_c = 299792458.0;         // Lichtgeschwindigkeit c
double const c_gamma = 1.0/(1+1/1.0);   // Der Exponent des Polytropen EOS mit dem Polytropenindex 0.5  rho = (P/K)^gamma    
double const c_k = 1.0;                 // Konstante K
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
    // Dies sind Funktionen, die einen Fehler hatten.
    //- (c_g(m + 4* c_pi*pow(r,3)*p/(pow(c_c,2)))*((pow((p/c_k),c_gamma))+p/pow(c,2)))/(r(r-2*c_g*m/pow(c,2)))
    // (-1)*(c_g*(m+4*c_pi*pow(r,3)*p/pow(c_c,2))*(pow((p/c_k),c_gamma))+p/(pow(c_c,2)))/((r*(r-c_g)*m)/pow(c_c,2))
    
    // hier wird die Funktion dm/dr definiert, welche für (r,p) einen Wert ausgibt
    double func2(double r, double p){
    double F2 = 4*c_pi*pow(p/c_k,c_gamma)*pow(r,2);
    return F2;
    }


int main(){
     
     
    ofstream file1;
    file1.open ("datenpunkte.csv"); //wird jetzt noch als csv abgespeichert, da noch keine Zeit gehabt um Gnuplot anzuschauen
        
        // In der While-Schleife wird das Runge-Kutta-Verfahren durchlaufen, bis der Druck 0 ist.
        while (pn > 1){
        K1f = c_delr*func1(rn, pn, mn); //Die K_i werden hier berechnet
        K1g = c_delr*func2(rn,pn);
        K2f = c_delr*func1(rn+ (1.0/2) *c_delr,pn+(1.0/2)*K1f,mn+(1.0/2)*K1g);
        K2g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K1f);
        K3f = c_delr*func1(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f,mn+(1.0/2)*K2g);
        K3g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f);
        K4f = c_delr*func1(rn+c_delr,pn+K3f,mn+K3g);
        K4g = c_delr*func2(rn+c_delr,pn+K3f);
        
        pn += (1.0/6) *K1f+(1.0/3)*K2f+(1.0/3)*K3f+(1.0/6)*K4f; //Der Wert p(r+delr) ergibt sich hier
        mn += (1.0/6) *K1g+(1.0/3)*K2g+(1.0/3)*K3g+(1.0/6)*K4g; //Der Wert m(r+delr) ergibt sich hier
        rn +=c_delr;
        // cout<<K1f<<endl;
        // cout<<rn<<"     "<<pn<<"   "<<mn<<endl; //hier werden die Zahlen überprüft, um den Fehler zu finden

     
        file1 <<rn<<","<<pn<<","<<mn<<endl;
      
    
    }
    file1.close();
    return 0;
}
