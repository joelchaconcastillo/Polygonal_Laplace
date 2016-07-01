/**
    Autor: Joel Chacón Castillo
    Fecha: 01/07/2016
    Descripcion: Este fichero contiene la implementacion de los metodos para
                 calcular las funciones de forma washpress y Laplace en c++.
**/
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
using namespace std;
/**
    Calcular las coordenadas del espacio canonico de cada nodo en
    el vecindario definido.
**/
void CoordenadasNodosPoligono(vector< vector<double> > &Coordenadas,int N)
{
    //pi =  acos(-1); de forma alternativa....
    double theta = (2.0*M_PI)/N;

    if(N == 4)
    {
        Coordenadas[0][0] = -1;
        Coordenadas[1][0] = -1;
        Coordenadas[0][1] = 1;
        Coordenadas[1][1] = -1;
        Coordenadas[0][2] = 1;
        Coordenadas[1][2] = 1;
        Coordenadas[1][3] = -1;
        Coordenadas[2][3] = 1;
        return;
    }

    for(int i = 0; i < N ;i++)
    {
        double Angulo = i*theta;
        Coordenadas[0][i] = cos(Angulo);
        Coordenadas[1][i] = sin(Angulo);
    }
}
/**
    Se realiza el calculo del circuncentro dados tres puntos
**/
void CalcularCircuncentro(double a1, double a2, double b1, double b2, double x, double y, double &v1, double &v2, double &v1x, double &v1y, double &v2x, double &v2y)
{
    double ToleranciaDeterminante = 1e-18;
    double D = (a1 - x)*(b2 - y) - (b1 - x)*(a2 - y);
    if(fabs(D) < ToleranciaDeterminante)
    {
        cout << "Tres puntos del triángulo son colineales";
        exit(0);
    }
    double Termino1 = 0.5*(((a1 - x)*(a1 + x) + (a2 - y)*(a2 + y))*(b2 - y) - ((b1 - x)*(b1 + x) + (b2 - y)*(b2 + y))*(a2 - y));
    v1 = Termino1/D;
    double Termino2 = 0.5*(((b1 - x)*(b1 + x) + (b2 - y)*(b2 + y))*(a1 - x) - ((a1 - x)*(a1 + x) + (a2 - y)*(a2 + y))*(b1 - x));
    v2 = Termino2/D;
    double Dx = a2 - b2;
    double Dy = b1 - a1;
    double Termino3 = 0.5*((b1 + a1)*(b1 - a1) + (b2 + a2)*(b2 - a2));
      v1x = (x - v1)*Dx/D;
      v1y = (Termino3 + y*Dx - v1*Dy)/D;
      v2x = (-Termino3 + x*Dy - v2*Dx)/D;
      v2y = (y - v2)*Dy/D;
}
/**
    Imprime en pantalla el valor de cada funcion de forma de laplace y su respectiva derivada,
    dada una coordenada en el espacio canonico y el numero de lado del poligono
**/
void InformacionFuncionFormaLaplace(double xi1, double xi2, int NLados)
{
    /**
        Generar el vector de coordenadas del espacio natural
    **/
    vector< vector<double> > Coordenadas(2, vector<double> (NLados,0) );
    vector< vector<double> > Phi(3, vector<double> (NLados,0) );
    vector<double> SumaDerivada(2);
    vector<double>CirCen(2,0);
    vector<double>Derivada(2,0);
    double Sumatoria = 0;
    CoordenadasNodosPoligono(Coordenadas, NLados);
    for(int i = 0; i < NLados; i++)
    {
        double a1 = Coordenadas[0][i];
        double a2 = Coordenadas[1][i];
        int IndexNodePrev = abs( (i+NLados-1)  ) % NLados ;
        int IndexNodeNext = (i+1)% NLados;
        double b1 = Coordenadas[0][IndexNodePrev];
        double b2 = Coordenadas[1][IndexNodePrev];
        double p1, p2,p1_xi1, p1_xi2, p2_xi1, p2_xi2,q1, q2, q1_xi1, q1_xi2, q2_xi1, q2_xi2;
        CalcularCircuncentro(b1, b2, a1, a2, xi1, xi2, p1, p2, p1_xi1, p1_xi2, p2_xi1, p2_xi2);

        b1 = Coordenadas[0][IndexNodeNext];
        b2 = Coordenadas[1][IndexNodeNext];
        CalcularCircuncentro(a1, a2, b1, b2, xi1, xi2, q1, q2, q1_xi1, q1_xi2, q2_xi1, q2_xi2);

        double sI = sqrt((p1 - q1)*(p1 - q1) + (p2 - q2)*(p2 - q2));

        double sI_xi1 = ((p1 - q1)*(p1_xi1 - q1_xi1) + (p2 - q2)*(p2_xi1 - q2_xi1))/sI;
        double sI_xi2 = ((p1 - q1)*(p1_xi2 - q1_xi2) + (p2 - q2)*(p2_xi2 - q2_xi2))/sI;
        double hI = sqrt((xi1 - a1)*(xi1 - a1) + (xi2 - a2)*(xi2 - a2));
        double hI_xi1 = (xi1 - a1)/hI;
        double hI_xi2 = (xi2 - a2)/hI;

         double sI_hI = sI/hI;


         Derivada[0] = (sI_xi1 - sI_hI*hI_xi1)/hI;
         Derivada[1] = (sI_xi2 - sI_hI*hI_xi2)/hI;

         Sumatoria+=sI/hI;
         SumaDerivada[0] += Derivada[0];
         SumaDerivada[1] += Derivada[1];
         Phi[0][i] = sI_hI;
         Phi[1][i] = Derivada[0];
         Phi[2][i] = Derivada[1];
    }

    for(int i = 0; i < NLados; i++)
    {
        Phi[0][i] = Phi[0][i] / Sumatoria;
        Phi[1][i] = (Phi[1][i]-Phi[0][i]*SumaDerivada[0]) / Sumatoria;
        Phi[2][i] = (Phi[2][i]-Phi[0][i]*SumaDerivada[1]) / Sumatoria;
    }

    for(int i = 0; i < NLados; i++)
    {
        cout << " Ni(x): " << Phi[0][i];
        cout << " dNi_x1(x): " << Phi[1][i];
        cout << " dNi_x2(x): " << Phi[2][i];
        cout << endl;
    }

}

/**
    Realiza el calulo del area del triangulo cuyos nodos son a,b,c
**/
double AreaTriangulo(double a1, double a2, double b1, double b2, double c1, double c2)
{
    return ((a1-c1)*(b2-c2) - (b1-c1)*(a1-c2)  )/2.0;
}
/**
    Realiza el calculo de la derivada del triangulo.
**/
double DerivadaAreaTriangulo(double a1, double a2, double b1, double b2, double c1, double c2)
{
    double Termino1 = ( (a1-c1)*b2 + (b2-c2)*a1 ) /2.0;
    double Termino2 = ( (b1-c1)*a2 + (a2-c2)*b1 ) /2.0;
    return Termino1-Termino2;
}
/**
    Realiza el calculo de la funcion de forma de Wachpress
    dada una coordenada y el numero de lados del poligono.
**/
void InformacionFuncionFormaWachspress(double xi1, double xi2, int NLados)
{
    /**
        Generar el vector de coordenadas del espacio natural
    **/
    vector< vector<double> > Coordenadas(2, vector<double> (NLados,0) );
    vector< vector<double> > Phi(3, vector<double> (NLados,0) );
    vector<double> SumaDerivada(2);
    vector<double>CirCen(2,0);
    vector<double>Derivada(2,0);
    double SumatoriaPesos = 0;
    CoordenadasNodosPoligono(Coordenadas, NLados);
    for(int i = 0; i < NLados; i++)
    {


        int IndexNodePrev = abs( (i+NLados-1)  ) % NLados ;
        int IndexNodeNext = (i+1)% NLados;
        double a1 = Coordenadas[0][IndexNodePrev];
        double a2 = Coordenadas[1][IndexNodePrev];

        double b1 = Coordenadas[0][i];
        double b2 = Coordenadas[1][i];

        double c1 = Coordenadas[0][IndexNodeNext];
        double c2 = Coordenadas[1][IndexNodeNext];

        double Area1 = AreaTriangulo(xi1, xi2,a1, a2, b1, b2 );
        double Area2 = AreaTriangulo(xi1, xi2,b1,b2,c1,c2);

        double Peso = 2*(1- (xi1*xi1) - ( xi2*xi2)  ) * pow(sin(M_PI/NLados),3)*cos(M_PI/NLados);
        Peso/=(Area1*Area2);
        SumatoriaPesos+=Peso;
        Derivada[0] = (-4*xi1*xi1)*pow(sin(M_PI/NLados),3)*cos(M_PI/NLados);
        Derivada[0] /=(Area1*Area2);
        Derivada[1] =  (-4*xi2*xi2)*pow(sin(M_PI/NLados),3)*cos(M_PI/NLados) ;
        Derivada[1] /=(Area1*Area2);
        SumaDerivada[0] += Derivada[0];
        SumaDerivada[1] += Derivada[1];

        Phi[0][i] = Peso;

        Phi[1][i] = Derivada[0];
        Phi[2][i] = Derivada[1];
    }

    for(int i = 0; i < NLados; i++)
    {
        Phi[0][i] = Phi[0][i] / SumatoriaPesos;
        Phi[1][i] = (Phi[1][i]-Phi[0][i]*SumaDerivada[0]) / SumatoriaPesos;
        Phi[2][i] = (Phi[2][i]-Phi[0][i]*SumaDerivada[1]) / SumatoriaPesos;
    }

    for(int i = 0; i < NLados; i++)
    {
        cout << " Ni(x): " << Phi[0][i];
        cout << " dNi_x1(x): " << Phi[1][i];
        cout << " dNi_x2(x): " << Phi[2][i];
        cout << endl;
    }

}
