#include <iostream>
#include <algorithm>
#include <random>
#include <math.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <ctime>
#include <functional>
#include <numeric>
#include <string>
#include <vector>
#include <climits>
#include "aco.h"
#include "graph.h"
#include "helpers.h"
#include <unistd.h>
using namespace std;


ACO::~ACO(){
    for (auto &hormiga : hormigas){
        for (Camion* c : hormiga.vector_camiones)
            delete c;
        hormiga.vector_camiones.clear();
    }
}


/* Constructor de la clase ACO
    Esta clase implementa la metaheuristica ACO para resolver el problema de recoleccion de residuos domiciliarios.
    El flujo de la metaheuristica es el siguiente:
    - Se inicializan las hormigas en un nodo inicial aleatorio dentro del conjunto de nodos iniciales.
    - Se inicializan las feromonas segun el parametro Tau, que es 1 por defecto TODO: hay una formula para establecer el mejor valor. ref: aco-book.
    - Se construye una solución para cada hormiga.
        + Se elige el siguiente nodo a visitar para cada hormiga.
        + Se visita el nodo elegido.
    - Se evapora las feromonas.
    - Se actualizan las feromonas.
    - Se repite el proceso hasta que se cumpla el criterio de parada.

    Parámetros:
        - graph: Puntero al grafo que se utilizará para resolver el problema
        - num_hormigas: Número de hormigas que se utilizarán para resolver el problema
*/
ACO::ACO(Graph *instancia, ACOArgs parametros_base)
{
    set_parametros(parametros_base);
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm *now_tm = std::localtime(&now_c);

    instancia->metadatos.encabezado["METODO"] = metodo == 0 ? "AntSystem" : metodo == 1 ? "MaxMin"
                                                                                        : "ACS";
    nombre_metodo = instancia->metadatos.encabezado["METODO"];
    std::string prefijo = instancia->metadatos.encabezado["NOMBRE"] + "-" + instancia->metadatos.encabezado["METODO"] + "-";
    nombre_instancia_salida = instancia->metadatos.encabezado["NOMBRE"];

    // Construir la ruta al archivo en la carpeta "output"
    std::strftime(filename, 100, "%Y%m%d%H%M%S", now_tm); // Formato: -AAAAMMDDHHMMSS

    //directorio_salida = ("Output/" + prefijo + std::string(filename));
    directorio_salida = parametros_base.directorio_salida;
    std::filesystem::create_directories(directorio_salida);
    if (!parametros_base.rescol)
        prefijo_salida = parametros_base.prefijo_salida + "_" + std::string(filename);
    else
        prefijo_salida = parametros_base.prefijo_salida;
    //directorio_salida.string() + "/" + parametros_base.prefijo_salida + "_" + std::string(filename) + ".txt";
    //set_filename(nombre_archivo_salida);

    // Inicializa las hormigas.
    for (int i = 0; i < num_hormigas; i++)
    {
        // Crea una nueva hormiga
        Hormiga hormiga;
        hormiga.id = i;

        // Establece la posición inicial de la hormiga de manera aleatorio entre la cantidad de nodos iniciales permitidos.
        hormiga.nodo_actual = nodoInicialAleatorio(instancia);

        // Inicializa un mapa para que la hormiga lleve la cuenta de las pasadas por los arcos.
        for (auto &par : instancia->arcos)
        {
            hormiga.arcos_visitados_tour[par.second] = 0;
            hormiga.arcos_visitados_salida[par.second] = 0;
        }

        hormiga.arcos_no_visitados = hormiga.arcos_visitados_tour;

        // Añade la hormiga a la lista de hormigas
        hormigas.push_back(hormiga);
    }

    // Establece el grafo.
    grafo = instancia;

    LimiteDeMejoras = parametros_base.LimiteDeMejoras;

    // Establece los arreglos para el algoritmo de dijkstra
    //dijkstra(*grafo, dist, next);

    int n = grafo->nodos.size();
    std::unordered_map<int, int> nodeToIndex;
    std::unordered_map<int, int> indexToNode;
    int index = 0;

    // Crear el mapeo de identificadores de nodos a índices consecutivos
    for (const auto &nodoPair : grafo->nodos) {
        nodeToIndex[nodoPair.first] = index;
        indexToNode[index] = nodoPair.first;
        index++;
    }
    dist.resize(n, std::vector<double>(n, std::numeric_limits<double>::infinity()));
    next.resize(n, std::vector<int>(n, -1));

    // Inicializar distancias y siguiente nodo
    for (const auto &nodoPair : grafo->nodos) {
        int u = nodeToIndex[nodoPair.first];
        dist[u][u] = 0;
        for (const auto &arco : nodoPair.second.saliente) {
            int v = nodeToIndex[arco.destino->id];
            dist[u][v] = arco.costo_recorrido;
            next[u][v] = v;
        }
    }

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }
}

/*
    Nodo aleatorio
    Genera un numero aleatorio entre 1 y el size del vector de nodos iniciales y retorna el nodo en aquel posicion.
*/
Nodo* ACO::nodoInicialAleatorio (Graph* instancia){

    int nodo_inicial = static_cast<int>(generar_numero_aleatorio(1, instancia->metadatos.nodos_iniciales.size() - 1));
    return &instancia->metadatos.nodos_iniciales[nodo_inicial];
}

/*
    Itera el algoritmo
    Este método itera el algoritmo, moviendo todas las hormigas, evaporando las feromonas y actualizando las feromonas.
*/
void ACO::iterar()
{
    // Mueve todas las hormigas.
    for (auto &hormiga : hormigas)
    {
        construirSolucion(hormiga);
        limpiar_rastro();        
        
        /*
        if (timeout_flag){
            continue;        
        }*/
    }
}

/*
    Construye la solución para una hormiga
    Este método construye la solución para una hormiga, moviendola por el grafo hasta que se cumpla el criterio de parada.

    Parámetros:
    - hormiga: Referencia a la hormiga que se está moviendo
*/
void ACO::construirSolucion(Hormiga &hormiga)
{
    Nodo *actual = nullptr;
    Nodo *siguiente = nullptr;
    acumulador_tiempo = 0;
    int aux = 0;
    //bool flag_camion = false;

    
    while (!solucionCompleta(hormiga))
    {
        auto start = std::chrono::high_resolution_clock::now();
        timeout_flag = false;

        //Inicializacion para la primera iteracion
        if (hormiga.vector_camiones.empty()){
            Camion *camion = new Camion();
            hormiga.vector_camiones.push_back(camion);
        } 

        actual = hormiga.nodo_actual;
        
        if (debug)
        {
        std::cout << "Hormiga numero " << hormiga.id << " en el nodo " << actual->id << std::endl;
        }
        if(sin_nuevo*2 > LimiteDeMejoras){
            //cout << "Entra a la mejora" << endl;
            std::vector<Nodo*> ruta_nueva = findPath(hormiga, ArcosNoVisitadoObligatorios(hormiga, 0));
            ruta_nueva.push_back(ArcosNoVisitadoObligatorios(hormiga, 1));
            sin_nuevo = 0;
            for(size_t i = 0; i < ruta_nueva.size(); i++) {
                if(hormiga.nodo_actual->id != ruta_nueva[i]->id) {
                    visitar(hormiga, ruta_nueva[i], aux);
                }
            }
        }else
            {
            siguiente = eligeSiguiente(hormiga, aux);
            if (!siguiente){
                if(usaDijkstra){
                    buscarDijkstra(hormiga,aux);
                    hormiga.copia_camino_tour.clear();
                }
                Camion *camion = new Camion();
                hormiga.vector_camiones.push_back(camion);
                hormiga.nodo_actual = &grafo->metadatos.nodos_iniciales[int((generar_numero_aleatorio(0, grafo->metadatos.nodos_iniciales.size() - 1)))];
                hormiga.vector_camiones[aux]->camino_final.insert(hormiga.vector_camiones[aux]->camino_final.end(), hormiga.vector_camiones[aux]->camino_tour.begin(), hormiga.vector_camiones[aux]->camino_tour.end());
                hormiga.vector_camiones[aux]->camino_final.insert(hormiga.vector_camiones[aux]->camino_final.end(), hormiga.vector_camiones[aux]->camino_salida.begin(), hormiga.vector_camiones[aux]->camino_salida.end());
                hormiga.vector_camiones[aux]->longitud_camino_final = hormiga.vector_camiones[aux]->longitud_camino_tour + hormiga.vector_camiones[aux]->longitud_camino_salida;
                hormiga.saltosSalida = hormiga.vector_camiones[aux]->longitud_camino_final - hormiga.vector_camiones[aux]->longitud_camino_tour;
                aux++;
                continue;
                } else visitar(hormiga, siguiente, aux);
            }
        /*std::cout << "esta la ruta de la hormiga " << std::endl;
        for(int i = 0; hormiga.camino_tour.size() > i ; i++){
            std::cout << hormiga.camino_tour[i].origen->id << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;*/
        if (usar_oscilador == 1)
        {
            oscilador.oscilar();
        } else if (usar_oscilador == 2){
            oscilador.oscilar_caotico();
        } 
        //cout << "alfa: " << ACO::alfa << " beta: " << ACO::beta << " rho: " << ACO::rho << " tau: " << ACO::tau << endl;
        auto stop = std::chrono::high_resolution_clock::now();        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        acumulador_tiempo += duration.count();
    }
    
    for (auto *camion : hormiga.vector_camiones)
    {
        if (camion && camion->camino_final.empty())
        {
            camion->camino_final.insert(camion->camino_final.end(), camion->camino_tour.begin(), camion->camino_tour.end());
            camion->camino_final.insert(camion->camino_final.end(), camion->camino_salida.begin(), camion->camino_salida.end());
            camion->longitud_camino_final = camion->longitud_camino_tour + camion->longitud_camino_salida;
            hormiga.saltosSalida = camion->longitud_camino_final - camion->longitud_camino_tour;
        }
    }

    file << "Epoca: " << epoca_actual << ", Evaluacion: " << evaluaciones << ", Mejor costo: " << mejor_costo << endl;
    // hormiga.copia_camino_final.insert(hormiga.copia_camino_final.end(), hormiga.copia_camino_tour.begin(), hormiga.copia_camino_tour.end());
    // hormiga.copia_camino_final.insert(hormiga.copia_camino_final.end(), hormiga.copia_camino_salida.begin(), hormiga.copia_camino_salida.end());
    
    calcular_costo_camino_camion(hormiga);
    hormiga.longitud_final_camiones += total_camino_tour_camion(hormiga) + total_camino_salida_camion(hormiga);
    hormiga.saltosSalida = hormiga.longitud_final_camiones - total_camino_tour_camion(hormiga); //no se para que es
    evaluaciones++;
}

double ACO::total_camino_tour_camion(Hormiga &hormiga){

    double total_longitud = 0;

    for (auto &camion : hormiga.vector_camiones){
        total_longitud += camion->longitud_camino_tour;
    }
    return total_longitud;

}

double ACO::total_camino_salida_camion(Hormiga &hormiga){
    
    double total_longitud = 0;

    for (auto &camion : hormiga.vector_camiones){
        total_longitud += camion->longitud_camino_salida;
    }
    return total_longitud;
}

#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <limits>

std::vector<Nodo*> ACO::findPath(Hormiga &hormiga, Nodo* targetNode) {
    std::vector<Nodo*> path;
    if (hormiga.nodo_actual == nullptr || targetNode == nullptr) {
        return path;
    }

    std::unordered_map<Nodo*, Nodo*> came_from;
    std::unordered_map<Nodo*, double> cost_so_far;
    auto compare = [&cost_so_far](Nodo* left, Nodo* right) { return cost_so_far[left] > cost_so_far[right]; };
    std::priority_queue<Nodo*, std::vector<Nodo*>, decltype(compare)> frontier(compare);

    came_from[hormiga.nodo_actual] = nullptr;
    cost_so_far[hormiga.nodo_actual] = 0;
    frontier.push(hormiga.nodo_actual);

    while (!frontier.empty()) {
        Nodo* current = frontier.top();
        frontier.pop();

        if (current == targetNode) {
            break;
        }

        for (Arco &arco : current->saliente) {
            if (hormiga.arcos_visitados_final.find(&arco) != hormiga.arcos_visitados_final.end() &&
                hormiga.arcos_visitados_final[&arco] > 0) {
                continue;
            }

            Nodo* next = arco.destino;
            double new_cost = cost_so_far[current] + arco.costo_recorrido;

            if (cost_so_far.find(next) == cost_so_far.end() || new_cost < cost_so_far[next]) {
                cost_so_far[next] = new_cost;
                came_from[next] = current;
                frontier.push(next);
            }
        }
    }

    if (came_from.find(targetNode) == came_from.end()) {
        return path;  // No path found
    }

    for (Nodo* current = targetNode; current != nullptr; current = came_from[current]) {
        path.push_back(current);
    }

    std::reverse(path.begin(), path.end());
    return path;
}


/*
    Elige el siguiente nodo
    Este método elige el siguiente nodo a visitar para una hormiga segun una regla de transicion pseudo aleatoria proporcional.
    El primer for recorre los nodos adyacentes al nodo actual de la hormiga
    El segundo for recorre los arcos buscando si el arco es el que conecta el nodo actual con el nodo adyacente.
    Cuando se encuentra el arco, se calcula la probabilidad de transición, se normaliza y se elige el nodo siguiente.

    Parámetros:
    - hormiga: Referencia a la hormiga que se está moviendo

    Retorna:
    - Nodo: Puntero al nodo elegido
*/
Nodo *ACO::eligeSiguiente(Hormiga &hormiga, int &aux)
{
    double total = 0.0;
    double r = generar_numero_aleatorio(0, 1.00);
    double acumulado = 0.0;
    double cantidad = 0.0;
    double tau_eta = 0.0;
    int pasadas = 0;
    int obligatorios_actuales = 0;
    //bool flag_falla = false;
    bool pasado_una_vez_por_obligatorias = false;
    Nodo *nodo = nullptr;
    std::unordered_map<Arco *, double> probabilidad;
    Metadatos metadatos;

    if (hormiga.nodo_actual == nullptr) {
        throw std::runtime_error("Nodo actual de hormiga es nulo.");
    }

    if (!hormiga.vector_camiones[aux]) throw std::runtime_error("Camión actual es nulo.");

    for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
    {
        if(i.first->obligatoria){
            obligatorios_actuales++;
            if(i.first->veces_recorrida > 0){
                pasadas++;
            }
        }
    }
    if(sin_nuevo >= LimiteDeMejoras && !hormiga.vector_camiones[aux]->camino_tour.empty()){
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            //flag_falla = true;
            if (!hormiga.vector_camiones[aux]->camino_tour.empty()){
                // si no esta vacio, se empiezan a comprobar los casos
                if (i.first->bidireccional == true)
                { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                    
                    // si el destino del reciproco no es el nodo actual, lo considero   
                    if (hormiga.copia_camino_tour.back().id != i.first->arco_reciproco->id && i.first->destino->id != hormiga.nodo_actual->id)
                    {
                        Arco *arco = nullptr;
                        arco = i.first;
                        if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante){ 
                            tau_eta = 1;
                            probabilidad[arco] = tau_eta;
                            total += tau_eta;
                            //cout << "if 1" << endl;
                            continue;
                        } else {
                            if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante){
                               //cout << "else 1" << endl;
                                return nullptr;
                            }
                        }
                    }
                }
                else
                { // si no es bidireccional, se agrega a las probabilidades de paso
                    Arco *arco = nullptr;
                    arco = i.first;
                    if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "if 2" << endl;
                        tau_eta = 1;
                        probabilidad[arco] = tau_eta;
                        total += tau_eta;
                    }else {
                        if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "else 2" << endl;
                            return nullptr;
                        }
                    }
                }
            }
        }
    }

    else if(pasado_una_vez_por_obligatorias == false){
        // No se pasa 1 vez por cada arista obligatoria
        // Se elige segun regla de transicion que arista obligatoria pasar
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            if (!hormiga.vector_camiones[aux]->camino_tour.empty())
            {
                // si no esta vacio, se empiezan a comprobar los casos
                if (i.first->bidireccional == true)
                { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                    // si el destino del reciproco no es el nodo actual, lo considero   
                    if (hormiga.copia_camino_tour.back().id != i.first->arco_reciproco->id && i.first->destino->id != hormiga.nodo_actual->id)
                    {
                        //Si el camion [aux] de la hormiga tiene capacidad y tiempo suficiente para pasar al siguiente nodo
                        Arco *arco = nullptr;
                        arco = i.first;
                        if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "if 3" << endl;
                            if (arco->veces_recorrida <= valor_limitador)
                            {
                                cantidad = hormiga.feromonas_locales[arco].cantidad;
                                if (arco->obligatoria == true){
                                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                                }
                                probabilidad[arco] = tau_eta;
                                total += tau_eta;
                                if (debug)
                                    cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                            }else {
                                continue;
                            }
                        } else if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "else 3" << endl;
                            return nullptr;
                        }
                    }
                }
                else
                { // si no es bidireccional, se agrega a las probabilidades de paso
                    
                    Arco *arco = nullptr;
                    arco = i.first;
                    if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "if 4" << endl;
                        if (arco->veces_recorrida <= valor_limitador)
                        {
                            cantidad = hormiga.feromonas_locales[arco].cantidad;
                            if (arco->obligatoria == true){
                                tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                            }
                            total += tau_eta;
                            probabilidad[arco] = tau_eta;
                        } else {
                            continue;
                        }
                    }else if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "else 4" << endl;
                        return nullptr;
                    }
                }
            }
            else
            {
                // si el camino esta vacio, simplemente se agrega el primero que encuentre sin comprobaciones extra
                Arco *arco = nullptr;
                arco = i.first;
                cantidad = hormiga.feromonas_locales[arco].cantidad;
                if (arco->obligatoria == true){
                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                }
                else {
                    tau_eta = 1;
                }
                probabilidad[arco] = tau_eta;
                total += tau_eta;
                if (debug){
                    cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                }
            }
        }
    }

    else if(pasado_una_vez_por_obligatorias == true){
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            if (!hormiga.vector_camiones[aux]->camino_tour.empty())
            {
                // si no esta vacio, se empiezan a comprobar los casos
                if (i.first->bidireccional == true)
                { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                    // si el destino del reciproco no es el nodo actual, lo considero   
                    if (hormiga.copia_camino_tour.back().id != i.first->arco_reciproco->id && i.first->destino->id != hormiga.nodo_actual->id)
                    {
                        Arco *arco = nullptr;
                        arco = i.first;
                        if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "if 5" << endl;
                            if (arco->veces_recorrida <= valor_limitador)
                            {
                                cantidad = hormiga.feromonas_locales[arco].cantidad;
                                if (arco->obligatoria == true){
                                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                                }
                                else {
                                    tau_eta = 1;
                                }
                                probabilidad[arco] = tau_eta;
                                total += tau_eta;
                                if (debug)
                                    cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                                    
                            }else {
                                continue;
                            }
                        } else if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "else 5" << endl;
                            return nullptr;
                        }
                    }
                }
                else
                { // si no es bidireccional, se agrega a las probabilidades de paso
                    Arco *arco = nullptr;
                    arco = i.first;
                    if (arco->costo_recoleccion<= hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "if 6" << endl;
                        if (arco->veces_recorrida <= valor_limitador)
                        {

                            cantidad = hormiga.feromonas_locales[arco].cantidad;
                            if (arco->obligatoria == true){
                                tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                            }
                            else {
                                tau_eta = 1;
                            }
                            total += tau_eta;
                            probabilidad[arco] = tau_eta;
                        } else {
                            continue;
                        }
                    } else if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "else 6" << endl;
                        return nullptr;
                    }
                }
            }
            else
            {
                // si el camino esta vacio, simplemente se agrega el primero que encuentre sin comprobaciones extra
                Arco *arco = nullptr;
                arco = i.first;
                cantidad = hormiga.feromonas_locales[arco].cantidad;
                if (arco->obligatoria == true){
                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                }
                else {
                    tau_eta = 1;
                }
                probabilidad[arco] = tau_eta;
                total += tau_eta;
                if (debug)
                    cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
            }
        }
    }

    if (total == 0){ // si nos quedamos sin pasadas se elige cualquiera segun el metodo normal
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id]){
            // si no esta vacio, se empiezan a comprobar los casos
            if (i.first->bidireccional == true)
            { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                if (hormiga.copia_camino_tour.back().id != i.first->arco_reciproco->id && i.first->destino->id != hormiga.nodo_actual->id)
                {
                    Arco *arco = nullptr;
                    arco = i.first;
                    if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "if 7" << endl;
                        cantidad = hormiga.feromonas_locales[arco].cantidad;
                        if (arco->obligatoria == true){
                            tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                        } 
                        else {
                            tau_eta = 1;
                        }
                        probabilidad[arco] = tau_eta;
                        total += tau_eta;
                        if (debug)
                            cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                    } else {
                        //cout << "else 7" << endl;
                        return nullptr;
                    }
                }  
            }
            else
            { // si no es bidireccional, se agrega a las probabilidades de paso
                Arco *arco = nullptr;
                arco = i.first;

                cantidad = hormiga.feromonas_locales[arco].cantidad;
                if (arco->obligatoria == true){
                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                }
                else {
                    tau_eta = 1;
                }
                tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                total += tau_eta;
                probabilidad[arco] = tau_eta;
            }
        }
    }
    if (total == 0){ // si nos quedamos sin pasadas se elige cualquiera segun el metodo normal
        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id]){
            // si no esta vacio, se empiezan a comprobar los casos
            
            if (i.first->bidireccional == true)
            { // si es bidireccional, se comprueba que el siguiente arco no sea una vuelta en U o si es la unica opcion, en este ultimo caso, se agrega a las probabilidades de paso
                Arco *arco = nullptr;
                arco = i.first;
                if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                    //cout << "if 8" << endl;
                    cantidad = hormiga.feromonas_locales[arco].cantidad;
                    if (arco->obligatoria == true){
                        tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                    } 
                    else {
                        tau_eta = 1;
                    }
                    probabilidad[arco] = tau_eta;
                    total += tau_eta;
                    if (debug)
                        cout << "arco:" << arco->origen->id << " " << arco->destino->id << " tau_eta: " << tau_eta << endl;
                    } else {
                        if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                            //cout << "else 8" << endl;
                            return nullptr;
                        }
                    }
            }
            else
            { // si no es bidireccional, se agrega a las probabilidades de paso
                Arco *arco = nullptr;
                arco = i.first;
                if (arco->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante ){
                    //cout << "if 9" << endl;
                    cantidad = hormiga.feromonas_locales[arco].cantidad;
                    if (arco->obligatoria == true){
                        tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                    }
                    else {
                        tau_eta = 1;
                    }
                    tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta);
                    total += tau_eta;
                    probabilidad[arco] = tau_eta;
                } else {
                    if (arco->costo_recoleccion > hormiga.vector_camiones[aux]->capacidad_restante ){
                        //cout << "else 9" << endl;
                        return nullptr;
                    }
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
            if(p.first->veces_recorrida == 0){
                sin_nuevo = 0;
            }
            else{
                sin_nuevo++;
            }

            break;
        }
    }

  
    //cout << "Arcos faltantes: " << hormiga.arcos_no_visitados.size() << " con r " << r << endl;
    if (debug)
        cout << "r: " << r << endl;
    if (debug)
        cout << "nodo elegido: " << nodo->id << endl;
    if (nodo == nullptr){
        /*if(flag_falla){
            cout << "Nodo actual sin salida: " << hormiga.nodo_actual->id<< endl;
            cout << "entrando a" << endl;
            cout << endl;
        }*/
        cout << "Nodo actual sin salida: " << hormiga.nodo_actual->id<< endl; // usar exepciones o algo try catch
        cout << "nodo nulo con r " << r << endl; // usar exepciones o algo try catch
    }
    //cout << "nodo elegido: " << nodo->id << endl;
    return nodo;

}


/*
    Visita el nodo siguiente
    Este método visita el nodo siguiente para una hormiga, actualizando el camino, el costo y la longitud del camino.

    Parámetros:
    - hormiga: Referencia a la hormiga que se está moviendo
    - nodo: Puntero al nodo que se va a visitar
*/
void ACO::visitar(Hormiga &hormiga, Nodo *nodo, int &aux)
{
    Arco *arco = nullptr;
    Camion camion;


    for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
    {
        if (i.first->destino->id == nodo->id)
        {
            arco = i.first;
            break;
        }
    }

    if (arco->veces_recorrida <= 1){
        hormiga.arcos_visitados_tour[arco] += 1;
    }
    
    hormiga.vector_camiones[aux]->camino_tour.push_back(*arco);
    hormiga.copia_camino_tour.push_back(*arco);
    hormiga.arcos_no_visitados.erase(arco);

    if (arco->bidireccional) {
        hormiga.arcos_no_visitados.erase(arco->arco_reciproco);
    }

    hormiga.nodo_actual = nodo;
    hormiga.vector_camiones[aux]->longitud_camino_tour += 1;
    
    if (arco->veces_recorrida == 0) {

        if (hormiga.feromonas_locales[arco].cantidad < umbral_inferior) hormiga.feromonas_locales[arco].cantidad = umbral_inferior;
        else hormiga.feromonas_locales[arco].cantidad *= (1 - rho_secundario);

        hormiga.vector_camiones[aux]->capacidad_restante -= (arco->costo_recoleccion);
        //hormiga.vector_camiones[aux]->tiempo_restante -= arco->costo_recorrido;
        hormiga.vector_camiones[aux]->costo_camino_camion += (arco->costo_recorrido + arco->costo_recoleccion);
    }
    else {

        if (hormiga.feromonas_locales[arco].cantidad < umbral_inferior) hormiga.feromonas_locales[arco].cantidad = umbral_inferior;
        else hormiga.feromonas_locales[arco].cantidad *= (1 - rho);
        //hormiga.vector_camiones[aux]->costo_camino_camion += arco->costo_recoleccion;
    }

    arco->veces_recorrida += 1;


    //cout << "El arco visitado es " << arco->origen->id << "-" << arco->destino->id<< " id: " << arco->id << endl;
    

    // if (i.first->costo_recoleccion <= hormiga.vector_camiones[aux]->capacidad_restante && i.first->destino->costo_recorrido <= hormiga.vector_camiones[aux].tiempo_restante) {
        
    //     hormiga.nodo_actual = nodo;
    // }

    
    /*double resta_feromona = rho;
    if(resta_feromona - 0.01*sin_nuevo >= 1){
        resta_feromona = 1;
    }

    double resta_feromona_secundario = rho_secundario;
    if(resta_feromona_secundario - 0.01*sin_nuevo >= 1){
        resta_feromona_secundario = 1;
    }*/
    if (arco->veces_recorrida == 1)
    {
        if (hormiga.feromonas_locales[arco].cantidad < umbral_inferior)
        {
            hormiga.feromonas_locales[arco].cantidad = umbral_inferior;
        }
        else
        {
            hormiga.feromonas_locales[arco].cantidad *= (1 - rho);
        }
        //hormiga.vector_camiones[aux]->costo_total += arco->costo_recoleccion;
    }
    else
    {
        if (hormiga.feromonas_locales[arco].cantidad < umbral_inferior)
        {
            hormiga.feromonas_locales[arco].cantidad = umbral_inferior;
        }
        else
        {
            hormiga.feromonas_locales[arco].cantidad *= (1 - rho_secundario);
        }
        //hormiga.vector_camiones[aux]->costo_total += arco->costo_recoleccion;
    }
    return;
}
// No se utiliza
void ACO::buscarSalidaMatriz(Hormiga &hormiga, int &aux)
{

    while (!enNodoTerminal(hormiga)) //
    {
        Nodo *nodo = nullptr;
        std::unordered_map<Arco *, double> probabilidad;
        double total = 0.0;
        double r = generar_numero_aleatorio(0, 1.00);
        double acumulado = 0.0;
        double cantidad = 0.0;
        double tau_eta = 0.0;

        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            Arco *arco = nullptr;
            arco = i.first;
            cantidad = feromonas_salida[arco].cantidad;
            tau_eta = pow(cantidad, alfa) * pow(grafo->informacion_heuristica[hormiga.nodo_actual->id][i.first], beta_salida);
            probabilidad[arco] = tau_eta;
            total += tau_eta;
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

        Arco *arco = nullptr;

        for (auto i : grafo->informacion_heuristica[hormiga.nodo_actual->id])
        {
            if (i.first->destino->id == nodo->id)
            {
                arco = i.first;
                break;
            }
        }

        arco->veces_recorrida += 1;
        hormiga.arcos_visitados_salida[arco] += 1;
        hormiga.vector_camiones[aux]->camino_salida.push_back(*arco);
        hormiga.nodo_actual = nodo;
        hormiga.vector_camiones[aux]->longitud_camino_salida += 1;
        //hormiga.vector_camiones[aux]->costo_camino_camion += (arco->costo_recoleccion/2);
    }
    return;
}

void ACO::buscarDijkstra(Hormiga &hormiga, int &aux)
{
    //std::vector<std::vector<double>> dist;
    //std::vector<std::vector<int>> next;
    //floydWarshall(*grafo, dist, next);

    std::unordered_map<int, int> nodeToIndex;
    std::unordered_map<int, int> indexToNode;
    int index = 0;

    // Crear el mapeo de identificadores de nodos a índices consecutivos
    for (const auto &nodoPair : grafo->nodos) {
        nodeToIndex[nodoPair.first] = index;
        indexToNode[index] = nodoPair.first;
        index++;
    }

    int nodoActualId = hormiga.nodo_actual->id;
    int nodoActualIndex = nodeToIndex[nodoActualId];
    int nodoSalidaIndex = -1;

    // Encontrar el nodo de salida más cercano
    double minDist = std::numeric_limits<double>::infinity();
    for (const auto &nodoSalida : grafo->metadatos.nodos_termino) {
        int terminoIndex = nodeToIndex[nodoSalida.id];
        if (dist[nodoActualIndex][terminoIndex] < minDist) {
            minDist = dist[nodoActualIndex][terminoIndex];
            nodoSalidaIndex = terminoIndex;
        }
    }

    if (nodoSalidaIndex == -1) {
        // No se encontró un nodo de salida válido
        return;
    }

    // Reconstruir la ruta más corta
    std::vector<int> path;
    int u = nodoActualIndex;
    while (u != nodoSalidaIndex) {
        u = next[u][nodoSalidaIndex];
        if (u == -1) {
            // No hay camino posible
            return;
        }
        path.push_back(u);
    }

    // Actualizar la hormiga con la ruta encontrada
    for (size_t i = 0; i < path.size(); ++i) {
        int vIndex = path[i];
        int v = indexToNode[vIndex];
        Arco *arco = nullptr;
        for (auto &a : hormiga.nodo_actual->saliente) {
            if (a.destino->id == v) {
                arco = &a;
                break;
            }
        }

        
        if (arco) {
            arco->veces_recorrida += 1;
            hormiga.arcos_visitados_salida[arco] += 1;
            hormiga.vector_camiones[aux]->camino_salida.push_back(*arco);
            hormiga.vector_camiones[aux]->longitud_camino_salida += 1;
            hormiga.nodo_actual = arco->destino;
            hormiga.vector_camiones[aux]->costo_camino_camion += (arco->costo_recorrido);
            if (enNodoTerminal(hormiga)) calcular_costo_camino_camion(hormiga);
        }
    }
    return;
}

void ACO::calcular_costo_camino_camion(Hormiga &hormiga){

    hormiga.costo_camino = 0.0;
    for (auto &camion : hormiga.vector_camiones){
        hormiga.costo_camino += camion->costo_camino_camion;
    }
}

void ACO::calcular_longitud_camiones(Hormiga &hormiga){

    for (auto &camion : hormiga.vector_camiones){
        if (!camion) continue;
        hormiga.longitud_final_camiones += camion->longitud_camino_tour + camion->longitud_camino_salida;
    }
}

int ACO::contarArcosNoVisitadosObligatorios(Hormiga &hormiga) {
    int count = 0;
    for (const auto &par : hormiga.arcos_no_visitados) {
        Arco *arco = par.first;
        if (arco->obligatoria) {
            count++;
        }
    }
    return count;
}

Nodo* ACO::ArcosNoVisitadoObligatorios(Hormiga &hormiga, int flag) {
    //int count = 0;
    Nodo* aux;
    for (const auto &par : hormiga.arcos_no_visitados) {
        Arco *arco = par.first;
        if (arco->obligatoria) {
            if(flag == 0){
                aux = arco -> origen;
                break;
            }
            else{
                aux = arco -> destino;
                break;
            }
       }
            //break;
        //}
    }
    return aux;
}


/*
    Verifica si la solución es completa
    Este método verifica si la solución es completa, es decir, si se han visitado todos los arcos al menos una vez y si el nodo actual es uno de los nodos finales.

    Parámetros:
    - hormiga: Referencia a la hormiga que se está moviendo

    Retorna:
    - bool: True si la solución es completa, False en caso contrario.
*/
bool ACO::solucionCompleta(Hormiga &hormiga)
{
    if(contarArcosNoVisitadosObligatorios(hormiga) == 0){
        //cout << '1' << endl;
        return true;
    }
    return false;
    

    
    /*bool completo = hormiga.arcos_no_visitados.empty();
    if (!usarMatrizSecundaria)
    {
        bool terminado = enNodoTerminal(hormiga);
        return (completo && terminado);
    }

    return completo;*/
}

/*
    Verifica si la hormiga esta en un nodo final

    Parámetros:
    - hormiga: Referencia a la hormiga que se está moviendo

    Retorna:
    - bool: True si la hormiga esta en un nodo final
*/
bool ACO::enNodoTerminal(Hormiga &hormiga)
{
    // El nodo actual es uno de los nodos finales
    bool termino = false;

    if(hormiga.nodo_actual == nullptr){
        return termino;
    }

    for (auto &nodo_final : grafo->metadatos.nodos_termino)
    {
        if (hormiga.nodo_actual->id == nodo_final.id)
        {
            termino = true;
            break;
        }
    }

    return termino;
}

/*
    Muestra la solución
    Este método muestra la solución, mostrando el camino por cada camion y el costo de la hormiga.
*/
void ACO::mostrar_solucion(bool show_solucion)
{
    for (int i = 0; i < 161; i++)
        cout << "-";
    cout << endl;
    if (show_solucion)
    {
        cout << "La solución es: " << endl;
        int index = 1;
        {
            for (auto &camion : mejor_solucion.vector_camiones){
                cout << "Camion " << index <<": ";

                // for (auto &arco : camion->camino_final)
                // {
                //     if (&arco != &camion->camino_final.back())
                //         cout << arco.origen->id << " -> ";
                //     //if (arco.origen->id == )
                if (!camion->camino_final.empty()){
                    cout << camion->camino_final.front().origen->id;

                    for (auto &arco : camion->camino_final){
                        cout << " -> " << arco.destino->id;
                    }
                }
                // cout << arco.destino->id << endl;
                cout << endl;
                index++;
            }
        }
    }
    cout << endl;
    cout << "✨🏆 Mejor resultado 🏆✨" << endl;
    cout << endl;
    cout << "Mejor hormiga: " << mejor_solucion.id << " 🐜🥇" << endl;
    if (timeout_flag_global)
    {
        cout << "Mejor costo: " << "inf" << " ⏩" << endl;
    } else {
        std::cout << "Mejor costo: " << std::fixed << std::setprecision(3) << mejor_solucion.costo_camino << " ⏩" << std::endl;
    }
    cout << "Mejor longitud: " << mejor_solucion.longitud_final_camiones << " ⚡" << endl;
}

// void ACO::exportar_solucion(std::chrono::microseconds duration, ACOArgs parametros_base)
// {
//     // suma de costos de recoleccion de los arcos
//     int suma_recoleccion = 0;
//     int suma_recorrer = 0;
//     int costo_pesos_pasada = 0;
//     for (auto &arco : grafo->arcos)
//     {
//         suma_recoleccion += arco.second->costo_recoleccion;
//     }
//     for (auto &arco : mejor_solucion.vector_camiones->camino_final)
//     {
//         suma_recorrer += arco.costo_recorrido;
//     }
//     costo_pesos_pasada = suma_recorrer - suma_recoleccion;

//     cout << "El directorio de salida:" << directorio_salida.string()<< endl;

//      //"resultados_csv.txt";
//     ruta_archivo_salida_csv = directorio_salida.string() + "/" + prefijo_salida + ".csv";
//     std::ofstream archivo_salida_csv(ruta_archivo_salida_csv);
//     archivo_salida_csv << nombre_instancia_salida << ",";
//     archivo_salida_csv << nombre_metodo << ",";
//     //archivo_salida_csv << costo_pesos_pasada << ",";
//     archivo_salida_csv << mejor_solucion.costo_camino << ",";
//     archivo_salida_csv << mejor_solucion.longitud_final_camiones << ",";
//     archivo_salida_csv << mejor_solucion.camino_tour.front().origen->id << ",";
//     if (mejor_solucion.camino_salida.empty())
//         archivo_salida_csv << mejor_solucion.camino_tour.back().destino->id << ",";
//     else
//         archivo_salida_csv << mejor_solucion.camino_salida.back().destino->id << ",";
//     archivo_salida_csv << duration.count() << ",";
//     archivo_salida_csv << mejor_solucion.longitud_camino_tour << ",";     //
//     archivo_salida_csv << mejor_solucion.saltosSalida << ","; //
//     archivo_salida_csv << grafo->arcos.size() << ",";
//     archivo_salida_csv << suma_recoleccion << ",";
//     archivo_salida_csv << suma_recorrer << endl;
//     archivo_salida_csv.close();

//     ruta_archivo_config_salida_csv = directorio_salida.string() + "/" + prefijo_salida + "_config.csv";
//     std::ofstream archivo_config_salida_csv(ruta_archivo_config_salida_csv);
//     archivo_config_salida_csv << usarMatrizSecundaria << ",";
//     archivo_config_salida_csv << usaDijkstra << ",";
//     archivo_config_salida_csv << parametros_base.oscilador << ",";
//     archivo_config_salida_csv << parametros_base.limitador << ",";
//     archivo_config_salida_csv << parametros_base.valor_limitador << ",";
//     archivo_config_salida_csv << parametros_base.beta_0 << ",";
//     archivo_config_salida_csv << alfa << ",";
//     archivo_config_salida_csv << beta << ",";
//     archivo_config_salida_csv << beta_salida << ",";
//     archivo_config_salida_csv << rho << ",";
//     archivo_config_salida_csv << rho_secundario << ",";
//     archivo_config_salida_csv << rho_salida << ",";
//     archivo_config_salida_csv << tau << ",";
//     if (usar_iteraciones){
//         archivo_config_salida_csv << parametros_base.iteraciones_max << ",";
//         archivo_config_salida_csv << "-1" << ",";
//     } else {
//         archivo_config_salida_csv << "-1" << ",";
//         archivo_config_salida_csv << parametros_base.evaluaciones_maximas << ",";
//     }
//     archivo_config_salida_csv << parametros_base.umbral_inferior << ",";
//     archivo_config_salida_csv << parametros_base.num_hormigas << ",";
//     archivo_config_salida_csv << parametros_base.epocas << ",";
//     if (parametros_base.full_aleatorio){
//         archivo_config_salida_csv << "-1" << ",";
//     } else {
//         archivo_config_salida_csv << parametros_base.semilla << ",";
//     }
//     archivo_config_salida_csv << parametros_base.umbral_superior << ",";
//     archivo_config_salida_csv << parametros_base.umbral_sin_mejora_limite << ",";
//     archivo_config_salida_csv << parametros_base.a << ",";
//     archivo_config_salida_csv << parametros_base.q_0 << ",";
//     archivo_config_salida_csv << parametros_base.csi << ",";
//     archivo_config_salida_csv << parametros_base.usar_iteraciones << ",";
//     archivo_config_salida_csv << parametros_base.usar_evaluaciones << ",";
//     archivo_config_salida_csv << parametros_base.irace << ",";
//     archivo_config_salida_csv << parametros_base.silence << ",";
//     archivo_config_salida_csv << parametros_base.full_aleatorio << endl;
//     archivo_config_salida_csv.close();
    
//     //if (parametros_base.irace || parametros_base.silence )
//     {
//         if (usar_bd){
//             std::stringstream ss;
//             ss << "python POSTDB.py " << ruta_archivo_salida_csv << " " << ruta_archivo_config_salida_csv;
//             std::string comando3 = ss.str();
//             std::system(comando3.c_str()); 
//         }
//     }
//     //"resultados.txt";
//     std::string ruta_archivo_salida = directorio_salida.string() + "/" + prefijo_salida + ".txt";
//     std::ofstream archivo_salida(ruta_archivo_salida);
//     //archivo_salida << "Mejor hormiga: " << mejor_solucion.id << endl;
//     archivo_salida << "Mejor longitud: " << mejor_solucion.longitud_final_camiones << endl;
//     archivo_salida << "Longitud tour: " << mejor_solucion.saltosTour << endl;     //
//     archivo_salida << "Longitud salida: " << mejor_solucion.saltosSalida << endl; //
//     archivo_salida << "Cantidad arcos: " << grafo->arcos.size() << endl;
//     archivo_salida << "Costo recoleccion: " << suma_recoleccion << endl;
//     archivo_salida << "Costo recorrer: " << suma_recorrer << endl;
//     archivo_salida << "Costo pesos pasada: " << costo_pesos_pasada << endl;
//     archivo_salida << "Mejor costo: " << mejor_solucion.costo_camino << endl;
//     archivo_salida << "Nodo inicio: " << mejor_solucion.camino_tour.front().origen->id << endl;
//     if (mejor_solucion.camino_salida.empty())
//         archivo_salida << "Nodo fin: " << mejor_solucion.camino_tour.back().destino->id << endl;
//     else
//         archivo_salida << "Nodo fin: " << mejor_solucion.camino_salida.back().destino->id << endl;
//     archivo_salida << "Tiempo de resolucion: " << duration.count() << " microsegundos" << endl;
//     archivo_salida << "La solución es:" << endl;
//     for (auto &arco : mejor_solucion.camino_final)
//     {
//         archivo_salida << arco.origen->id << endl;
//     }
//     archivo_salida << mejor_solucion.camino_final.back().destino->id << endl;
//     archivo_salida.close();

//     cout << "La solución es: " << endl;
//     // Abre el archivo para escribir
//     std::string camino_file = directorio_salida.string() + "/" + prefijo_salida  + "_camino.txt";
//     std::ofstream archivo(camino_file);
//     // Escribe cada arco y la cantidad de veces que se pasó por él
//     for (const auto &arco : mejor_solucion.camino_final)
//     {
//         archivo << arco.origen->id << endl;
//     }
//     archivo << mejor_solucion.camino_final.back().destino->id << endl;
//     archivo.close();
// }

void ACO::exportar_mapa_resultados()
{
    for (auto &par : mejor_solucion.arcos_visitados_tour)
    {
        Arco *arco = par.first;
        int suma = par.second;

        // Sumar el valor del segundo mapa
        if (mejor_solucion.arcos_visitados_salida.find(arco) != mejor_solucion.arcos_visitados_salida.end())
        {
            suma += mejor_solucion.arcos_visitados_salida[arco];
        }

        // Agregar al mapa de suma
        mejor_solucion.arcos_visitados_final[arco] = suma;
    }
    // Abre el archivo para escribir
    std::string mapa_file = directorio_salida.string() + "/" + prefijo_salida  + "_mapa.txt";
    std::ofstream archivo(mapa_file);
    // Escribe cada arco y la cantidad de veces que se pasó por él
    // sumar arcos visitados de salida
    for (const auto &arco : mejor_solucion.arcos_visitados_final)
    {
        auto it = mejor_solucion.arcos_visitados_final.find(arco.first);
        if (it != mejor_solucion.arcos_visitados_final.end())
        {
            archivo << "(" << arco.first->origen->id << ", " << arco.first->destino->id << ") "
                    << " : " << (arco.second) << endl;
        }
    }

    archivo.close();
}
/*
    Este método guarda la mejor solución, es decir, la hormiga con el mejor(menor) costo y longitud de camino.

    Retorna:
    - Hormiga: La mejor hormiga
*/

void ACO::exportar_solucion(std::chrono::microseconds duration, ACOArgs parametros_base)
{
    long long suma_recoleccion = 0;
    long long suma_recorrer = 0;
    long long costo_pesos_pasada = 0;
    double total_longitud_tour = 0.0;
    double total_longitud_salida = 0.0;
    for (auto &arco : grafo->arcos)
    {
        suma_recoleccion += arco.second->costo_recoleccion;
    }

    Nodo *nodo_inicio = nullptr;
    Nodo *nodo_fin = nullptr;

    for (auto *camion : mejor_solucion.vector_camiones)
    {
        if (!camion)
            continue;

        total_longitud_tour += camion->longitud_camino_tour;
        total_longitud_salida += camion->longitud_camino_salida;

        for (auto &arco : camion->camino_final)
        {
            suma_recorrer += arco.costo_recoleccion;
        }

        if (!camion->camino_final.empty())
        {
            if (!nodo_inicio)
            {
                nodo_inicio = camion->camino_final.front().origen;
            }
            nodo_fin = camion->camino_final.back().destino;
        }
    }

    costo_pesos_pasada = suma_recorrer - suma_recoleccion;

    ruta_archivo_salida_csv = directorio_salida.string() + "/" + prefijo_salida + ".csv";
    std::ofstream archivo_salida_csv(ruta_archivo_salida_csv);
    archivo_salida_csv << nombre_instancia_salida << ",";
    archivo_salida_csv << nombre_metodo << ",";

    if (timeout_flag_global)
    {
        archivo_salida_csv << "inf" << ",";
    }
    else
    {
        std::ostringstream costo;
        costo << std::fixed << std::setprecision(3) << mejor_solucion.costo_camino;
        archivo_salida_csv << costo.str() << ",";
    }

    archivo_salida_csv << mejor_solucion.longitud_final_camiones << ",";
    archivo_salida_csv << (nodo_inicio ? nodo_inicio->id : -1) << ",";
    archivo_salida_csv << (nodo_fin ? nodo_fin->id : -1) << ",";
    archivo_salida_csv << duration.count() << ",";
    archivo_salida_csv << total_longitud_tour << ",";
    archivo_salida_csv << total_longitud_salida << ",";
    archivo_salida_csv << grafo->arcos.size() << ",";
    archivo_salida_csv << suma_recoleccion << ",";
    archivo_salida_csv << suma_recorrer << std::endl;
    archivo_salida_csv.close();

    ruta_archivo_config_salida_csv = directorio_salida.string() + "/" + prefijo_salida + "_config.csv";
    std::ofstream archivo_config_salida_csv(ruta_archivo_config_salida_csv);
    archivo_config_salida_csv << usarMatrizSecundaria << ",";
    archivo_config_salida_csv << usaDijkstra << ",";
    archivo_config_salida_csv << parametros_base.oscilador << ",";
    archivo_config_salida_csv << parametros_base.limitador << ",";
    archivo_config_salida_csv << parametros_base.valor_limitador << ",";
    archivo_config_salida_csv << parametros_base.beta_0 << ",";
    archivo_config_salida_csv << alfa << ",";
    archivo_config_salida_csv << beta << ",";
    archivo_config_salida_csv << beta_salida << ",";
    archivo_config_salida_csv << rho << ",";
    archivo_config_salida_csv << rho_secundario << ",";
    archivo_config_salida_csv << rho_salida << ",";
    archivo_config_salida_csv << tau << ",";
    if (usar_iteraciones)
    {
        archivo_config_salida_csv << parametros_base.iteraciones_max << ",";
        archivo_config_salida_csv << "-1" << ",";
    }
    else
    {
        archivo_config_salida_csv << "-1" << ",";
        archivo_config_salida_csv << parametros_base.evaluaciones_maximas << ",";
    }
    archivo_config_salida_csv << parametros_base.umbral_inferior << ",";
    archivo_config_salida_csv << parametros_base.num_hormigas << ",";
    archivo_config_salida_csv << parametros_base.epocas << ",";
    if (parametros_base.full_aleatorio)
    {
        archivo_config_salida_csv << "-1" << ",";
    }
    else
    {
        archivo_config_salida_csv << parametros_base.semilla << ",";
    }
    archivo_config_salida_csv << parametros_base.umbral_superior << ",";
    archivo_config_salida_csv << parametros_base.umbral_sin_mejora_limite << ",";
    archivo_config_salida_csv << parametros_base.a << ",";
    archivo_config_salida_csv << parametros_base.q_0 << ",";
    archivo_config_salida_csv << parametros_base.csi << ",";
    archivo_config_salida_csv << parametros_base.usar_iteraciones << ",";
    archivo_config_salida_csv << parametros_base.usar_evaluaciones << ",";
    archivo_config_salida_csv << parametros_base.irace << ",";
    archivo_config_salida_csv << parametros_base.silence << ",";
    archivo_config_salida_csv << parametros_base.full_aleatorio << std::endl;
    archivo_config_salida_csv.close();

    if (usar_bd)
    {
        std::stringstream ss;
        ss << "python POSTDB.py " << ruta_archivo_salida_csv << " " << ruta_archivo_config_salida_csv;
        std::string comando3 = ss.str();
        std::system(comando3.c_str());
    }

    std::string ruta_archivo_salida = directorio_salida.string() + "/" + prefijo_salida + ".txt";
    std::ofstream archivo_salida(ruta_archivo_salida);
    archivo_salida << "Mejor longitud: " << mejor_solucion.longitud_final_camiones << std::endl;
    archivo_salida << "Longitud tour total: " << total_longitud_tour << std::endl;
    archivo_salida << "Longitud salida total: " << total_longitud_salida << std::endl;
    archivo_salida << "Cantidad arcos: " << grafo->arcos.size() << std::endl;
    archivo_salida << "Costo recoleccion: " << suma_recoleccion << std::endl;
    archivo_salida << "Costo recorrer: " << suma_recorrer << std::endl;
    archivo_salida << "Costo pesos pasada: " << costo_pesos_pasada << std::endl;
    if (timeout_flag_global)
    {
        archivo_salida << "Mejor costo: inf" << std::endl;
    }
    else
    {
        std::ostringstream costo;
        costo << std::fixed << std::setprecision(3) << mejor_solucion.costo_camino;
        archivo_salida << "Mejor costo: " << costo.str() << std::endl;
    }
    archivo_salida << "Nodo inicio: " << (nodo_inicio ? nodo_inicio->id : -1) << std::endl;
    archivo_salida << "Nodo fin: " << (nodo_fin ? nodo_fin->id : -1) << std::endl;
    archivo_salida << "Tiempo de resolucion: " << duration.count() << " microsegundos" << std::endl;

    if (!mejor_solucion.vector_camiones.empty())
    {
        archivo_salida << "La solución por camión es:" << std::endl;
        int indice_camion = 1;
        for (auto *camion : mejor_solucion.vector_camiones)
        {
            if (!camion)
                continue;

            archivo_salida << "Camion " << indice_camion++ << ": ";
            if (!camion->camino_final.empty())
            {
                archivo_salida << camion->camino_final.front().origen->id;
                for (auto &arco : camion->camino_final)
                {
                    archivo_salida << " -> " << arco.destino->id;
                }
            }
            archivo_salida << std::endl;
            archivo_salida << "  Longitud tour: " << camion->longitud_camino_tour << std::endl;
            archivo_salida << "  Longitud salida: " << camion->longitud_camino_salida << std::endl;
            archivo_salida << "  Longitud total: " << camion->longitud_camino_final << std::endl;
            std::ostringstream costo_camion;
            costo_camion << std::fixed << std::setprecision(3) << camion->costo_camino_camion;
            archivo_salida << "  Costo recorrido: " << costo_camion.str() << std::endl;
        }
    }
    else
    {
        archivo_salida << "No se registraron camiones en la mejor solución." << std::endl;
    }
    archivo_salida.close();

    std::string camino_file = directorio_salida.string() + "/" + prefijo_salida + "_camino.txt";
    std::ofstream archivo_camino(camino_file);
    int indice_camion = 1;
    for (auto *camion : mejor_solucion.vector_camiones)
    {
        if (!camion)
            continue;

        archivo_camino << "Camion " << indice_camion++ << ":";
        if (!camion->camino_final.empty())
        {
            archivo_camino << " " << camion->camino_final.front().origen->id;
            for (auto &arco : camion->camino_final)
            {
                archivo_camino << " -> " << arco.destino->id;
            }
        }
        archivo_camino << std::endl;
    }
    archivo_camino.close();
}

Hormiga ACO::guardar_mejor_solucion_iteracion()
{
    /*
        - Construir funcion para sumar todos los costos_total de cada camion por hormiga.
        - modificar esta loop agregando costos de los camiones.
    */
    for (auto &hormiga : hormigas)
    {
        calcular_costo_camino_camion(hormiga);
        if (hormiga.solucion_valida){
            if (hormiga.costo_camino < mejor_costo && (hormiga.longitud_final_camiones) < mejor_longitud)
            {
                mejor_costo = hormiga.costo_camino;
                mejor_longitud = hormiga.longitud_final_camiones;
                mejor_hormiga = hormiga;
            }
        } else {
            continue;
        }
    }
    return mejor_hormiga;
}

/*
    Limpia la memoria y datos del algoritmo
    Este método limpia la memoria y datos del algoritmo, dejando el algoritmo en su estado inicial.
*/
void ACO::limpiar()
{
    for (auto &hormiga : hormigas)
    {

        hormiga.copia_camino_tour.clear();
        hormiga.vector_camiones.clear();
        // hormiga.camino_salida.clear();
        // hormiga.camino_final.clear();
        hormiga.solucion_valida = true;
        for (auto &par : grafo->arcos)
        {
            hormiga.arcos_visitados_tour[par.second] = 0;
            hormiga.arcos_visitados_salida[par.second] = 0;
        }
        hormiga.arcos_no_visitados = hormiga.arcos_visitados_tour;
        hormiga.longitud_final_camiones = 0;
        // hormiga.longitud_camino_tour = 0;
        // hormiga.longitud_camino_salida = 0;
        // hormiga.longitud_final_camiones = 0;
        hormiga.saltosSalida = 0;
        hormiga.saltosTour = 0;
        // limpiar vector de camiones.
        hormiga.costo_camino = 0;
        hormiga.feromonas_locales = feromonas;
        hormiga.nodo_actual = &grafo->metadatos.nodos_iniciales[int((generar_numero_aleatorio(0, grafo->metadatos.nodos_iniciales.size() - 1)))];
    }
}
void ACO::limpiar_rastro()
{
    for (auto &arco : grafo->arcos)
        arco.second->veces_recorrida = 0;
}

void ACO::reset()
{
    iteraciones = 0;
    limpiar();
    for (auto &arco : grafo->arcos)
        arco.second->veces_recorrida = 0;
}

// void ACO::limpiar_camiones(Hormiga &hormiga)
// {
//     for (auto *camion : hormiga.vector_camiones) 
//         delete camion;
//     hormiga.vector_camiones.clear();
// }

/*void ACO::abrir_file()
{
    file.open(nombre_archivo_salida);
}*/

/*void ACO::cerrar_file()
{
    file.close();
}*/

std::string ACO::get_filename()
{
    std::string nombre_archivo_salida = directorio_salida.string() + "/" + prefijo_salida + ".txt";
    return nombre_archivo_salida;
}

/*void ACO::set_filename(std::string filename)
{
    this->nombre_archivo_salida = filename;
}*/

void ACO::set_parametros(const ACOArgs parametros_base)
{
    nombre_instancia = parametros_base.nombre_instancia;
    metodo = parametros_base.metodo;

    alfa = parametros_base.alfa;
    float alfa_inc = parametros_base.inc_alfa;
    float alfa_min = parametros_base.min_alfa;
    float alfa_max = parametros_base.max_alfa;
    oscilador.agregarParametro(alfa,alfa_inc,alfa_min,alfa_max);

    
    beta_salida = parametros_base.beta_salida;
    //oscilador.agregarParametro(beta_salida,0.01,1,5); // el oscilador rompe la matriz de salida, no tocar


    rho = parametros_base.rho;
    float rho_inc = parametros_base.inc_rho;
    float rho_min = parametros_base.min_rho;
    float rho_max = parametros_base.max_rho;
    oscilador.agregarParametro(rho,rho_inc,rho_min,rho_max);

    rho_secundario = parametros_base.rho_secundario;
    float rho_secundario_inc = parametros_base.inc_rho_secundario;
    float rho_secundario_min = parametros_base.min_rho_secundario;
    float rho_secundario_max = parametros_base.max_rho_secundario;
    oscilador.agregarParametro(rho_secundario,rho_secundario_inc,rho_secundario_min,rho_secundario_max); 

    rho_salida = parametros_base.rho_salida;
    //oscilador.agregarParametro(rho_salida,0.01,0.1,0.5); // el oscilador rompe la matriz de salida, no tocar

    tau = parametros_base.tau;
    //oscilador.agregarParametro(tau,0.01,1,5);
    umbral_inferior = parametros_base.umbral_inferior;
    num_hormigas = parametros_base.num_hormigas;
    epocas = parametros_base.epocas;
    usarMatrizSecundaria = parametros_base.usaMatrizSecundaria;
    usaDijkstra = parametros_base.usaDijkstra;
    epoca_actual = 0;
    usar_iteraciones = parametros_base.usar_iteraciones;
    usar_evaluaciones = parametros_base.usar_evaluaciones;
    usar_tiempo = parametros_base.usar_tiempo;
    if ((usar_iteraciones && usar_evaluaciones) ||
        (usar_iteraciones && usar_tiempo) ||
        (usar_evaluaciones && usar_tiempo)) {
        std::cout << "No se puede usar más de un criterio de parada al mismo tiempo" << std::endl;
        exit(1);
    } else{
        if (usar_iteraciones){
            iteraciones_max = parametros_base.iteraciones_max;        
        } else if (usar_evaluaciones){
            evaluaciones_maximas = parametros_base.evaluaciones_maximas;
            iteraciones_max = evaluaciones_maximas / num_hormigas;
        } else if (usar_tiempo){
            comienzo_ejecucion = std::chrono::high_resolution_clock::now();
            tiempo_total_ejecucion = std::chrono::duration<double>(parametros_base.tiempo_maximo);
        } else{
            cout << "No se ha especificado un criterio de parada" << endl;
            //exit(1);
            //Para debug
            cout << "Usando iteraciones" << endl;
            iteraciones_max = parametros_base.iteraciones_max;        
            cout << iteraciones_max << endl;

        }
    } 
    usar_oscilador = parametros_base.oscilador;
    usar_bd = parametros_base.conectar_bd;
    beta = parametros_base.beta;
    bool beta0 = parametros_base.beta_0;
    if (beta0) {
        beta = 0;     
    } else {
        float beta_inc = parametros_base.inc_beta;
        float beta_min = parametros_base.min_beta;
        float beta_max = parametros_base.max_beta;
        oscilador.agregarParametro(beta,beta_inc,beta_min,beta_max); 
    }
    bool usar_limitador = parametros_base.limitador;
    if (usar_limitador){
        valor_limitador = parametros_base.valor_limitador;        
    } else {
        valor_limitador = INT_MAX;
    }


    
}

Hormiga ACO::get_mejor_solucion()
{
    return mejor_solucion;
}


