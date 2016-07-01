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
#include "FuncionesForma.hpp"
using namespace std;

int main()
{
    cout << "Funcion de Laplace:" << endl;
    InformacionFuncionFormaLaplace(0.3, 0.3, 7);
    cout << "Funcion de Wachspress:"<<endl;
    InformacionFuncionFormaWachspress(0.3, 0.3, 7);
    return 0;
}
