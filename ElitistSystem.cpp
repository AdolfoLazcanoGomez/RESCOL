
#include "ElitistSystem.h"
#include "aco.h"
#include <iostream>
#include <chrono>
#include "helpers.h"


ElitistSystem::ElitistSystem(Graph *instancia, ACOArgs parametros) : ACO(instancia, parametros)
{
    set_parametrosACS(parametros);
    inicializar_feromonas();
}

/* Resuelve el problema
    Este método resuelve el problema simplemente iterando el algoritmo hasta que se cumpla el criterio de parada.
    Algunas alternativas de mejora son:
    - Establecer un criterio de parada basado en la calidad de las soluciones, mientras menos mejor.
    - Establecer un criterio de parada basado en el tiempo de ejecución.
    - Establecer un criterio de parada basado en la cantidad de iteraciones sin mejora, de manera local, global y/o por hormiga.
    + Mas opciones en aco-book.
*/
void ElitistSystem::resolver()
{
    while (iteraciones < iteraciones_max)
    {
        iterar();
        mejor_solucion = guardar_mejor_solucion_iteracion();
        limpiar();
        iteraciones++;
    }
};

void ElitistSystem::iterar()
{
    ACO::iterar();
    // evaporar feromonas
    for (auto i = feromonas.begin(); i != feromonas.end(); i++)
    {
        if (i->second.cantidad < umbral_inferior)
        {
            i->second.cantidad = umbral_inferior;
        }
        else
        {
            i->second.cantidad *= (1 - rho);
        }        
    }
    if (usarMatrizSecundaria){
        for (auto i = feromonas_salida.begin(); i != feromonas_salida.end(); i++)
        {
            if (i->second.cantidad < umbral_inferior)
            {
                i->second.cantidad = umbral_inferior;
            }
            else
            {
                i->second.cantidad *= (1 - rho_salida);
            }        
        }
    }
    mejor_solucion = guardar_mejor_solucion_iteracion();
    // Actualiza las feromonas.

    if (mejor_solucion.solucion_valida){
        for (auto &par : mejor_solucion.arcos_visitados_tour)
        {
            Arco *a = par.first;
            int pasadas = par.second;
            if (pasadas != 0)
                feromonas.at(a).cantidad += (tau / (rho * pow(mejor_solucion.longitud_final_camiones, 2) * pasadas));
        }
        // if (usarMatrizSecundaria)
        // {
        //     for (auto &par : mejor_solucion.arcos_visitados_salida)
        //     {
        //         Arco *a = par.first;
        //         int pasadas = par.second;
        //         if (pasadas != 0)
        //             feromonas_salida.at(a).cantidad += (tau / (rho * pow(mejor_solucion.longitud_camino_salida, 2) * pasadas));
        //     }
        // }
    }
    
}





/*
    Guarda la mejor solución de la iteración
    Este método guarda la mejor solución de la iteración, para luego ser utilizada en la actualización de feromonas.
*/
void ElitistSystem::inicializar_feromonas()
{

    for (auto &par : grafo->arcos)
    {
        Arco *arco = par.second;
        Feromona feromona_inicial = {arco->origen, arco->destino, tau};
        feromonas_salida[arco] = feromona_inicial;
        feromonas[arco] = feromona_inicial;
        for (auto &hormiga : hormigas)
            hormiga.feromonas_locales[arco] = feromona_inicial;
    }
}

void ElitistSystem::set_parametrosACS(ACOArgs parametros)
{
    //q_0 = parametros.q_0;
    //csi = parametros.csi;
    alfa = 1;
    tau = 1; // grafo->nodos.size(); // TODO: Revisar esto
}

