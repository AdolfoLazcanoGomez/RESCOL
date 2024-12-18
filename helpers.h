#ifndef HELPERS_H
#define HELPERS_H
#include <string>
#include <random>
#include "aco.h"
using namespace std;

extern std::mt19937 generador;
void inicializar_generador(int semilla);
std::string eliminarEspacios(std::string str);
double generar_numero_aleatorio(double min = 1, double max = 999999999.00);
bool camino_tour_vacio(const Hormiga &hormiga);


#endif
