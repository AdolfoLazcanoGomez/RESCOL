#include <random>
#include "graph.h"
#include "aco.h"
#include <algorithm>


std::mt19937 generador;
void inicializar_generador(int semilla) {
    generador.seed(semilla);
}

// Función para eliminar los espacios de una cadena
std::string eliminarEspacios(std::string str) {
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());      
    return str;
}

// Esta función genera un número aleatorio entre min y max (ambos incluidos).
double generar_numero_aleatorio(double min, double max) {
    std::uniform_real_distribution<> distribucion(min, max);
    return distribucion(generador);
}

bool camino_tour_vacio(const Hormiga &hormiga){
    for (const auto& camion : hormiga.vector_camiones){
        if (!camion->camino_tour.empty()){
            return false;
        }
    }
    return true;
}

