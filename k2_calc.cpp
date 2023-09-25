/* Dieses Script rechnet die Punkte für die Relation M(R) für gegebenes P_c aus.*/
//test


#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

// dp/dr
double func1(double r, double p, double m); // hier werden die Funktionen deklariert
// dm/dr
double func2(double r, double p); 
//dH/dr
double func3(double b);
//db/dr
double func4(double p, double m, double r, double b, double h);

void funcmr(double p0);

double const c_g = 6.6743*pow(10,-11);       // Gravitationskonstante G
//double pn = 1.164*pow(10,31);              // Dies ist der Anfangswert des Drucks P(0)
double mn = 0.0;                             // Der Anfangswert der eingeschlossenen Masse m(0)
double rn = pow(10,-30);                     // Der Startwert des Radius
double bn = 2*rn;                            // Der Startwert für b
double hn = pow(rn,2);                       // Der Startwert für H
double const c_delr = 1.;                     // Die Schrittweite delta r
double const c_c = 299792458.0;              // Lichtgeschwindigkeit c
double nindex=1.;
double const c_gamma = (1+1.0/nindex);          // Der Exponent des Polytropen EOS mit dem Polytropenindex n, rho = (P/K)^gamma    
double const c_k = 3.56109*pow(10,8);             
double const c_pi= 3.14159265359;            // pi
double solmass = 1.98847e30; //Sonnenmasse in SI
double solmass_geom = solmass * pow(c_c, -2) * c_g;






// hier werden die Koeffizienten des Runge-Kutta-Verfahrens deklariert
double K1f; // f steht hier für die DGL dp/dr
double K1g; // g steht hier für die DGL dm/dr 
double K1h; // h steht hier für die DGL dh/dr
double K1b; // b steht hier für die DGL db/dr
double K2f;
double K2g;
double K2h;
double K2b;
double K3f;
double K3g;
double K3h;
double K3b;
double K4f;
double K4g;
double K4h;
double K4b;
double y; // r*b(R)/H(R)
double k2; // love number k2. Hieraus kann die Tidal Deformability berechnet werden
double c; //  M/R= Kompaktheit des Sterns





    // hier wird die Funktion dp/dr definiert, welche für (r,p,m) einen Wert ausgibt.
    double func1(double r, double p, double m){
        double F1 = - ((pow((p/c_k),1/c_gamma))+p)*(m+4*c_pi*pow(r,3)*p)/(r*(r-2*m));
        return F1;
    }
    // hier wird die Funktion dm/dr definiert, welche für (r,p) einen Wert ausgibt.
    double func2(double r, double p){
        double F2 = 4*c_pi*pow(r,2)*pow((p/c_k),1/c_gamma);
        return F2;
    }
    double M, R, bet, H, M2;

    //dH/dr
    double func3(double b){
        double F3 = b;
        return F3; 
    }
    
    //db/dr
    double func4(double p, double m, double r, double b, double h){
        double F4 = -b *(2/r+pow(1-2*m/r,-1)*(2*m/pow(r,2)+4*c_pi*r*(p-pow(p/c_k,1/c_gamma))))
                    -h*(-6*pow(1-2*m/r,-1)/pow(r,2)+4*c_pi*pow(1-2*m/r,-1)*(5*pow(p/c_k,1/c_gamma)
                    +9*p+(p+pow(p/c_k,1/c_gamma)/(c_gamma*c_k*pow(p/c_k,1-1/c_gamma))))
                    -4*pow((m+4*c_pi*pow(r,3)*p)/(r*(r-2*m)),2));//Hinderer 2008
        
        
        /*double F4 = 2*pow(1-2*m/r,-1)*h*(-2*c_pi*(5*pow(p/c_k,1/c_gamma)+9*p
        +(p+pow(p/c_k,1/c_gamma)/(c_gamma*c_k*pow(p/c_k,1-1/c_gamma))))+3/pow(r,2)
        +2*pow(1-2*m/r,-1)*pow(m/pow(r,2)+4*c_pi*r*p,2))
        +2*b/r*pow(1-2*m/r,-1)*(-1+m/r+2*c_pi*pow(r,2)*(pow(p/c_k,1/c_gamma)-p));*/ //Hinderer2010       
        return F4;
    }

//Diese Funktion weist M und R ihre Werte für ein gegebens P_c zu.
void funcmr(double pn){

    rn=pow(10,-30);
    mn=0.0;
    hn=pow(rn,2);
    bn=2*rn;
    while (pn > 1*pow(10,-18)){
        
        K1f = c_delr*func1(rn, pn, mn); //Die K_i werden hier berechnet (Kif->dp/dr Kig->dm/dr Kih->dH/dr Kib->db/dr)
        K1g = c_delr*func2(rn,pn);
        K1h = c_delr*func3(bn);
        K1b = c_delr*func4(pn, mn, rn, bn, hn);
        K2f = c_delr*func1(rn+ (1.0/2) *c_delr,pn+(1.0/2)*K1f,mn+(1.0/2)*K1g);
        K2g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K1f);
        K2h = c_delr*func3(bn+1.0/2*K1b);
        K2b = c_delr*func4(pn+1.0/2*K1f,mn+1.0/2*K1g,rn+1.0/2*c_delr,bn+1.0/2*K1b,hn+1.0/2*K1h);
        K3f = c_delr*func1(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f,mn+(1.0/2)*K2g);
        K3g = c_delr*func2(rn+(1.0/2)*c_delr,pn+(1.0/2)*K2f);
        K3h = c_delr*func3(bn+1.0/2*K2b);
        K3b = c_delr*func4(pn+1.0/2*K2f,mn+1.0/2*K2g,rn+1.0/2*c_delr,bn+1.0/2*K2b,hn+1.0/2*K2h);
        K4f = c_delr*func1(rn+c_delr,pn+K3f,mn+K3g);
        K4g = c_delr*func2(rn+c_delr,pn+K3f);
        K4h = c_delr*func3(bn+K3b);
        K4b = c_delr*func4(pn+K3f,mn+K3g,rn+c_delr,bn+K3b,hn+K3h);
        
        M= mn;
        R= rn;
        bet = bn;
        H = hn; 
        pn += (1.0/6) *K1f+(1.0/3)*K2f+(1.0/3)*K3f+(1.0/6)*K4f; //Der Wert p(r+delr) ergibt sich hier
        mn += (1.0/6) *K1g+(1.0/3)*K2g+(1.0/3)*K3g+(1.0/6)*K4g; //Der Wert m(r+delr) ergibt sich hier
        hn += (1.0/6) *K1h+(1.0/3)*K2h+(1.0/3)*K3h+(1.0/6)*K4h; //Der Wert h(r+delr) ergibt sich hier
        bn += (1.0/6) *K1b+(1.0/3)*K2b+(1.0/3)*K3b+(1.0/6)*K4b; //Der Wert b(r+delr) ergibt sich hier
        rn +=c_delr;

        
      
    
    }
}


// In der Main fct wird funcmr() aufgerufen und die Datenpunkte werde abgespeichert.
int main(){
     ofstream file1;
     M2=0;
    file1.open ("tidal_n1_st.dat");
     for (double p = 1.164*pow(10,31) *c_g/pow(c_c,4) ; p < 7.562*pow(10,45)*c_g/pow(c_c,4); p = p*1.1)
     {
        funcmr(p);
        if (/*M != 0 &&*/ M>M2)
        {
            
        
        
        y = R*bet/H;
        c = M/R;


        k2 =  (8*pow(c,5)/5)*pow(1-2*c,2)*(2+2*c*(y-1)-y)*
        pow(2*c*(6-3*y+3*c*(5*y-8))+
        4*pow(c,3)*(13-11*y+c*(3*y-2)+2*pow(c,2)*(1+y))+
        3*pow(1-2*c,2)*(2-y+2*c*(y-1))*log(1-2*c),-1);


        file1 <<M/solmass_geom<<'\t'<<2./3*k2*pow(c,-5)<<endl;
        M2=M;
        }
        else
        {
            p=7.562*pow(10,46)*c_g/pow(c_c,4);
        }
     }
     file1.close();
    return 0;
}
