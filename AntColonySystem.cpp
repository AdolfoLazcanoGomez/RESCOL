
#include "AntColonySystem.h"
#include "aco.h"
#include <iostream>
#include <chrono>
#include "helpers.h"


AntColonySystem::AntColonySystem(Graph *instancia, ACOArgs parametros) : ACO(instancia, parametros)
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
void AntColonySystem::resolver()
{
    if(usar_tiempo){
        while (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - comienzo_ejecucion) <= tiempo_total_ejecucion )
        {
            std::cout << "hola" << std::endl;
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
    /*int i = 0;
    while (i < 100)
    {
        i++;
        iterar();
        mejor_solucion = guardar_mejor_solucion_iteracion();
        limpiar();
        iteraciones++;
    }*/
};

void AntColonySystem::iterar()
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
        if (usarMatrizSecundaria)
        {
            for (auto &par : mejor_solucion.arcos_visitados_salida)
            {
                Arco *a = par.first;
                int pasadas = par.second;
                if (pasadas != 0)
                    feromonas_salida.at(a).cantidad += (tau / (rho * pow(mejor_solucion.longitud_final_camiones, 2) * pasadas));
            }
        }
    }
    
}


Nodo *AntColonySystem::eligeSiguiente(Hormiga &hormiga, int &aux)
{
    double total = 0.0;
    double r = generar_numero_aleatorio(0, 1.00);
    double q = generar_numero_aleatorio(0, 1.00);
    double acumulado = 0.0;
    double cantidad = 0.0;
    double cantidad2 = 0.0;
    double tau_eta = 0.0;

    Nodo *nodo = nullptr;
    std::unordered_map<Arco *, double> probabilidad;
    if (q < q_0)
    {
        double cantidad_max = 0;
        // Seleccionar el arco con mas feromonas 
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id]){
            cantidad = i.second;
            cantidad2 = feromonas[i.first].cantidad;
            double calculo = cantidad2 * pow(cantidad,beta);

            if (cantidad2 > cantidad_max){
                cantidad_max = cantidad;
                nodo = i.first->destino;
            }
        }
    }
    else
    {
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            if (!hormiga.vector_camiones[aux]->camino_tour.empty())
            {
                // si no esta vacio, se empiezan a comprobar los casos
                if (i.first->bidireccional == true)
                { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                    if ((hormiga.copia_camino_tour.back().id != i.first->arco_reciproco->id) || (grafo->informacion_heuristica[hormiga.nodo_actual->id].size() == 1))
                    {
                        if (i.first->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante && i.first->costo_recorrido <= hormiga.vector_camiones[aux]->tiempo_restante){
                            Arco *arco = nullptr;
                            arco = i.first;
                            // if (arco->veces_recorrida <= 4)
                            {
                                cantidad = hormiga.feromonas_locales[arco].cantidad;
                                tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                                probabilidad[arco] = tau_eta;
                                total += tau_eta;
                                if (debug)
                                    cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                            }
                        } else {
                            if (i.first->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante || i.first->costo_recorrido > hormiga.vector_camiones[aux]->tiempo_restante){
                                break; //return nodo
                            }
                        }
                    }
                }
                else
                { // si no es bidireccional, se agrega a las probabilidades de paso
                    if (i.first->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante && i.first->costo_recorrido <= hormiga.vector_camiones[aux]->tiempo_restante){
                        Arco *arco = nullptr;
                        arco = i.first;
                        // if (arco->veces_recorrida <= 4)
                        {

                            cantidad = hormiga.feromonas_locales[arco].cantidad;
                            tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                            total += tau_eta;
                            probabilidad[arco] = tau_eta;
                        }
                    } else {
                        if (i.first->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante || i.first->costo_recorrido > hormiga.vector_camiones[aux]->tiempo_restante){
                            break; //return nodo
                        }
                    }
                }
            }
            else
            {
                // si el camino esta vacio, simplemente se agrega el primero que encuentre sin comprobaciones extra
                if (i.first->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante && i.first->costo_recorrido <= hormiga.vector_camiones[aux]->tiempo_restante){
                    Arco *arco = nullptr;
                    arco = i.first;
                    cantidad = hormiga.feromonas_locales[arco].cantidad;
                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                    probabilidad[arco] = tau_eta;
                    total += tau_eta;
                    if (debug){
                        cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                    }
                } else {
                    if (i.first->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante || i.first->costo_recorrido > hormiga.vector_camiones[aux]->tiempo_restante){
                        break; //return nodo
                    }
                }
            }
        }
        for (auto &p : probabilidad)
        {
            acumulado += p.second / total;
            if (r <= acumulado)
            {
                nodo = p.first->destino;
                break;
            }
        }
        if (debug)
            cout << "r: " << r << endl;
        if (debug)
            cout << "nodo elegido: " << nodo->id << endl;
    }

    for (auto& pair : hormiga.feromonas_locales)
    {
        pair.second.cantidad = (1 - csi) * pair.second.cantidad + csi * tau;
    }

    if (nodo == nullptr)
    {
        cout << "No hay nodos disponibles" << endl;
    }
    return nodo;
}


/*
    Guarda la mejor solución de la iteración
    Este método guarda la mejor solución de la iteración, para luego ser utilizada en la actualización de feromonas.
*/
void AntColonySystem::inicializar_feromonas()
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

void AntColonySystem::set_parametrosACS(ACOArgs parametros)
{
    q_0 = parametros.q_0;
    csi = parametros.csi;
    alfa = parametros.alfa;
    tau = parametros.tau; // grafo->nodos.size(); // TODO: Revisar esto
}

