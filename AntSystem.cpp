#include "AntSystem.h"
#include "aco.h"
#include <chrono>
#include <atomic>
#include <thread>


AntSystem::AntSystem(Graph *instancia, ACOArgs parametros) : ACO(instancia, parametros)
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
void AntSystem::resolver()
{
    if(usar_tiempo){
        std::atomic<bool> time_up(false);
        while (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - comienzo_ejecucion) <= tiempo_total_ejecucion )
        {

            std::thread iterar_thread([&time_up, this]() {
                this->iterar();
                time_up = true; // Marca que la función ha terminado
            });

            std::thread timer_thread([&time_up, this]() {
                auto tiempo_total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(this->tiempo_total_ejecucion).count();
                for (int i = 0; i < tiempo_total_ms / 100 && !time_up && solucionCompleta(mejor_solucion) == true; ++i) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    //reemplazar_hormiga();
                }
                time_up = true;
            });

            iterar_thread.join();
            timer_thread.join();
            
            //iterar();
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
void AntSystem::iterar()
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
    // if (usarMatrizSecundaria){
    //     for (auto i = feromonas_salida.begin(); i != feromonas_salida.end(); i++)
    //     {
    //         if (i->second.cantidad < umbral_inferior)
    //         {
    //             i->second.cantidad = umbral_inferior;
    //         }
    //         else
    //         {
    //             i->second.cantidad *= (1 - rho_salida);
    //         }        
    //     }
    // }
    // Actualiza las feromonas.
    for (Hormiga &hormiga : hormigas)
    {
        if (hormiga.solucion_valida){
            for (auto &par : hormiga.arcos_visitados_tour)
            {
                Arco *a = par.first;
                int pasadas = par.second;
                if (pasadas != 0)
                    feromonas.at(a).cantidad += (tau / (rho * pow(hormiga.longitud_final_camiones, 2) * pasadas));
                    //cout << endl;
            }
            // if (usarMatrizSecundaria)
            // {
            //     for (auto &par : hormiga.arcos_visitados_salida)
            //     {
            //         Arco *a = par.first;
            //         int pasadas = par.second;
            //         if (pasadas != 0)
            //             feromonas_salida.at(a).cantidad += (tau / (rho * pow(hormiga.longitud_final_camiones, 2) * pasadas));
            //     }
            // }
        }
    }
}

/*
    Guarda la mejor solución de la iteración
    Este método guarda la mejor solución de la iteración, para luego ser utilizada en la actualización de feromonas.
*/
void AntSystem::inicializar_feromonas()
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

void AntSystem::reemplazar_hormiga(){ 
    for(int i = 0; i < hormigas.size(); i++){ 
        if(solucionCompleta(hormigas[i])){ 
            mejor_solucion = hormigas[i];
            return; 
        } 
    } 
}