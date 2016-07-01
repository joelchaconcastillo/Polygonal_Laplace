#ifndef FUNCIONESFORMA_HPP_INCLUDED
#define FUNCIONESFORMA_HPP_INCLUDED
using namespace std;
/**
    Calcular las coordenadas del espacio canonico de cada nodo en
    el vecindario definido.
**/
void CoordenadasNodosPoligono(vector< vector<double> > &Coordenadas,int N);
/**
    Se realiza el calculo del circuncentro dados tres puntos
**/
void CalcularCircuncentro(double a1, double a2, double b1, double b2, double x, double y, double &v1, double &v2, double &v1x, double &v1y, double &v2x, double &v2y);
/**
    Imprime en pantalla el valor de cada funcion de forma de laplace y su respectiva derivada,
    dada una coordenada en el espacio canonico y el numero de lado del poligono
**/
void InformacionFuncionFormaLaplace(double xi1, double xi2, int NLados);
/**
    Realiza el calulo del area del triangulo cuyos nodos son a,b,c
**/
double AreaTriangulo(double a1, double a2, double b1, double b2, double c1, double c2);
/**
    Realiza el calculo de la derivada del triangulo.
**/
double DerivadaAreaTriangulo(double a1, double a2, double b1, double b2, double c1, double c2);
/**
    Realiza el calculo de la funcion de forma de Wachpress
    dada una coordenada y el numero de lados del poligono.
**/
void InformacionFuncionFormaWachspress(double xi1, double xi2, int NLados);


#endif // FUNCIONESFORMA_HPP_INCLUDED
