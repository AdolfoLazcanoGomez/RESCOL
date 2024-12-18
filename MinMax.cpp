#include "MinMax.h"
#include "aco.h"
#include <iostream>
#include <chrono>
#include <ctime>

MinMax::MinMax(Graph *instancia, ACOArgs parametros) : ACO(instancia, parametros)
{
    set_parametros(parametros);
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
void MinMax::resolver()
{
    int iteraciones_sin_mejora = 0;
    if(usar_tiempo){
        while (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - comienzo_ejecucion) <= tiempo_total_ejecucion)
        {        
            iterar();
            Hormiga mejor_solucion_anterior = mejor_solucion;
            mejor_solucion = guardar_mejor_solucion_iteracion();

            if (mejor_solucion_anterior.costo_camino == mejor_solucion.costo_camino && mejor_solucion_anterior.longitud_final_camiones == mejor_solucion.longitud_final_camiones)
            {
                iteraciones_sin_mejora++;
                if (iteraciones_sin_mejora == umbral_sin_mejora_limite)
                {
                    inicializar_feromonas();
                    iteraciones_sin_mejora = 0;
                }
            }
            else
            {
                umbral_superior = tau / (rho * mejor_solucion.longitud_final_camiones);
                umbral_inferior = umbral_superior / a;
                iteraciones_sin_mejora = 0;
            }
            limpiar();
            iteraciones++;
        }
    }
    else{
        while (iteraciones < iteraciones_max)
        {        
            iterar();
            Hormiga mejor_solucion_anterior = mejor_solucion;
            mejor_solucion = guardar_mejor_solucion_iteracion();

            if (mejor_solucion_anterior.costo_camino == mejor_solucion.costo_camino && mejor_solucion_anterior.longitud_final_camiones == mejor_solucion.longitud_final_camiones)
            {
                iteraciones_sin_mejora++;
                if (iteraciones_sin_mejora == umbral_sin_mejora_limite)
                {
                    inicializar_feromonas();
                    iteraciones_sin_mejora = 0;
                }
            }
            else
            {
                umbral_superior = tau / (rho * mejor_solucion.longitud_final_camiones);
                umbral_inferior = umbral_superior / a;
                iteraciones_sin_mejora = 0;
            }
            limpiar();
            iteraciones++;
        }
    }
};

void MinMax::iterar()
{
    ACO::iterar();

    // evaporar feromonas
    for (auto i = feromonas.begin(); i != feromonas.end(); i++)
    {
        i->second.cantidad *= (1 - rho);
        if (i->second.cantidad < umbral_inferior)
            i->second.cantidad = umbral_inferior;
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

    // Actualiza las feromonas.
    for (auto &par : mejor_solucion.arcos_visitados_tour)
    {
        Arco *a = par.first;
        int pasadas = par.second;
        if (pasadas != 0)
            feromonas.at(a).cantidad += (tau / (rho * mejor_solucion.longitud_final_camiones * pasadas));
        if (feromonas.at(a).cantidad > umbral_superior)
            feromonas.at(a).cantidad = umbral_superior;
    }
    if (usarMatrizSecundaria)
    {   
        for (auto &par : mejor_solucion.arcos_visitados_salida)
        {
            Arco *a = par.first;
            int pasadas = par.second;
            if (pasadas != 0)
                feromonas_salida.at(a).cantidad += (tau / (rho * mejor_solucion.longitud_final_camiones * pasadas));
            if (feromonas_salida.at(a).cantidad > umbral_superior)
                feromonas_salida.at(a).cantidad = umbral_superior;
        }
    }
}

/*
    Guarda la mejor solución de la iteración
    Este método guarda la mejor solución de la iteración, para luego ser utilizada en la actualización de feromonas.
*/
void MinMax::inicializar_feromonas()
{
    umbral_inferior = umbral_superior_inicial / a;
    umbral_superior = umbral_superior_inicial;
    for (auto &par : grafo->arcos)
    {
        Arco *arco = par.second;
        Feromona feromona_inicial = {arco->origen, arco->destino, umbral_superior_inicial};
        feromonas[arco] = feromona_inicial;
        feromonas_salida[arco] = feromona_inicial;
        for (auto &hormiga : hormigas)
            hormiga.feromonas_locales[arco] = feromona_inicial;
    }
}

void MinMax::set_parametros(ACOArgs parametros)
{
    umbral_superior = parametros.umbral_superior;
    umbral_superior_inicial = parametros.umbral_superior;
    umbral_sin_mejora_limite = parametros.umbral_sin_mejora_limite;
    a = parametros.a;
}
