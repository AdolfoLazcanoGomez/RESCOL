#include "ElitistSystem2.h"
#include "aco.h"
#include <chrono>

ElitistSystem2::ElitistSystem2(Graph *instancia, ACOArgs parametros) : ACO(instancia, parametros)
{
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
void ElitistSystem2::resolver()
{

    if(usar_tiempo){
        std::cout << "elitist" << std::endl;
        while (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - comienzo_ejecucion) <= tiempo_total_ejecucion )
        {
            //std::cout << "holi1" << std::endl;
            iterar();
            mejor_solucion = guardar_mejor_solucion_iteracion();
            limpiar();
            iteraciones++;
        }
    }
    else{
        while (iteraciones < iteraciones_max)
        {
            iterar();
            mejor_solucion = guardar_mejor_solucion_iteracion();
            limpiar();
            iteraciones++;
        }
    }
}
void ElitistSystem2::iterar()
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

    // Actualiza las feromonas.

    if (mejor_solucion.solucion_valida){
        for (auto &par : mejor_solucion.arcos_visitados_tour)
        {
            Arco *a = par.first;
            int pasadas = par.second;
            if (pasadas != 0)
                feromonas.at(a).cantidad += elitist*(tau / (rho * pow(mejor_solucion.longitud_final_camiones, 2) * pasadas));
        }
        if (usarMatrizSecundaria)
        {
            for (auto &par : mejor_solucion.arcos_visitados_salida)
            {
                Arco *a = par.first;
                int pasadas = par.second;
                if (pasadas != 0)
                    feromonas_salida.at(a).cantidad += elitist*(tau / (rho * pow(mejor_solucion.longitud_final_camiones, 2) * pasadas));
            }
        }
    }





    

    for (Hormiga &hormiga : hormigas)
    {
        if (hormiga.solucion_valida){
            for (auto &par : hormiga.arcos_visitados_tour)
            {
                Arco *a = par.first;
                int pasadas = par.second;
                if (pasadas != 0)
                    feromonas.at(a).cantidad += (tau / (rho * pow(hormiga.longitud_final_camiones, 2) * pasadas));
            }
            if (usarMatrizSecundaria)
            {
                for (auto &par : hormiga.arcos_visitados_salida)
                {
                    Arco *a = par.first;
                    int pasadas = par.second;
                    if (pasadas != 0)
                        feromonas_salida.at(a).cantidad += (tau / (rho * pow(hormiga.longitud_final_camiones, 2) * pasadas));
                }
            }
        }
    }
}

/*
    Guarda la mejor solución de la iteración
    Este método guarda la mejor solución de la iteración, para luego ser utilizada en la actualización de feromonas.
*/
void ElitistSystem2::inicializar_feromonas()
{
    for (auto &par : grafo->arcos)
    {
        Arco *arco = par.second;
        Feromona feromona_inicial = {arco->origen, arco->destino, tau};
        feromonas[arco] = feromona_inicial;
        feromonas_salida[arco] = feromona_inicial;
        for (auto &hormiga : hormigas)
            hormiga.feromonas_locales[arco] = feromona_inicial;
    }
}


void ElitistSystem2::set_parametrosACS(ACOArgs parametros)
{
    elitist = parametros.factor_elitista;
}