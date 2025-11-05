#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

#include <thread>
#include <chrono>

#include "argparse/argparse.hpp"

#include "graph.h"
#include "aco.h"
#include "AntSystem.h"
#include "MinMax.h"
#include "AntColonySystem.h"
#include "ElitistSystem.h"
#include "ElitistSystem2.h"
#include "instance_reader.h"
#include "helpers.h"
#include "enums.h"
#include "test.h"


#include <cstdio>
#include <iostream>
#include <string>

using namespace std;

/*
    Implementacion de la metaheuristica ACO para el problema de recoleccion de residuos.
    Autor: Alonso M. Robledo Álamos.
    Fecha: 2023-07-23
    Correo electrónico: alonsorobledo@protonmail.com | alonso.robledo.a@mail.pucv.cl

    Este programa implementa una metaheuristica tipo Ant Colony Optimization (ACO) para resolver un problema de ruteo de arcos (ARP)
    en el contexto de la recoleccion de residuos domiciliarios.
    Esto en el contexto de mi proyecto de tesis para optar al grado de magister en ingenieria informatica.
*/

// Funcion principal
/* Parámetros:
    - nombre_archivo: Nombre del archivo que contiene la instancia
    - iteraciones: Número de itecout << "Iteracion: " << iteraciones << endl;raciones que realizará el algoritmo
    - hormigas: Número de hormigas que se utilizarán
    - leer_restricciones: Indica si se leerán las restricciones del archivo
    - debug: Indica si se mostrarán mensajes de debug
*/
int main(int argc, char *argv[])
{


    auto comienzo = std::chrono::high_resolution_clock::now();

    ConfigPrograma config;
    Graph grafo = Graph();
    ACO *aco;
    ACOArgs parametros_base = argparse::parse<ACOArgs>(argc, argv);
    if (parametros_base.full_aleatorio){
        inicializar_generador(std::random_device{}());  
    } else {
        inicializar_generador(parametros_base.semilla);  
    }
    try {
        std::cout << "Instancia recibida: " << parametros_base.nombre_instancia << std::endl;
        grafo = leerInstancia(parametros_base.nombre_instancia, config.leer_restricciones, config.leer_coordenadas);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Fallo al leer la instancia: " << e.what() << std::endl;
        return 1;
    }

    if (config.debug)
    {
        cout << endl;
        for (int i = 0; i < 161; i++)
            cout << "-";
        cout << endl;
        cout << "Datos de debug 🪲: " << endl;
        cout << endl;
        cout << "Nombre de la instancia: " << grafo.metadatos.encabezado["NOMBRE"] << endl;
        cout << "Comentario: " << grafo.metadatos.encabezado["COMENTARIO"] << endl;
        cout << "Número de nodos: " << grafo.metadatos.encabezado["VERTICES"] << endl;
        cout << "Número de arcos obligatorios: " << grafo.metadatos.encabezado["ARISTAS_REQ"] << endl;
        (grafo.metadatos.encabezado["ARISTAS_NOREQ"] != "0") ? cout << "Número de arcos opcionales: " << grafo.metadatos.encabezado["ARISTAS_NOREQ"] << endl : cout << "No hay arcos opcionales" << endl;
        cout << "Número de nodos iniciales: " << grafo.metadatos.encabezado["NODOS_INICIALES"] << endl;
        (grafo.metadatos.encabezado["NODOS_TERMINO"] != "0") ? cout << "Número de nodos terminales: " << grafo.metadatos.encabezado["NODOS_TERMINO"] << endl : cout << "Todos los nodos pueden ser terminales" << endl;
        cout << "Número de nodos terminales: " << grafo.metadatos.encabezado["NODOS_TERMINO"] << endl;
        (grafo.metadatos.encabezado["RESTRICCIONES"] != "0") ? cout << "Número de restricciones: " << grafo.metadatos.encabezado["RESTRICCIONES"] << endl : cout << "No hay restricciones" << endl;
        cout << endl;
        cout << "Datos: " << endl;
        cout << "Nodos iniciales: " << endl;
        for (auto &nodo : grafo.metadatos.nodos_iniciales)
        {
            cout << nodo.id << " " << endl;
        }
        cout << endl;
        cout << "Nodos terminales: " << endl;
        for (auto &nodo : grafo.metadatos.nodos_termino)
        {
            cout << nodo.id << " ";
        }
        cout << endl;
        for (int i = 0; i < 161; i++)
            cout << "-";
        cout << endl;
    }
    cout << endl;
    switch (parametros_base.metodo)
    {
    case ANT_SYSTEM:
        aco = new AntSystem(&grafo, parametros_base);
        break;
    case MIN_MAX:
        aco = new MinMax(&grafo, parametros_base);
        break;
    case ELITIST:
        aco = new ElitistSystem2(&grafo, parametros_base);
        break;
    case ACS:
        aco = new AntColonySystem(&grafo, parametros_base);
        break;
    case TEST:
        bfs(grafo.metadatos.nodos_iniciales[0].id, grafo.informacion_heuristica);
    }
    ///////////////////////
    if (!grafo.metadatos.nodos_iniciales.empty()) {
    bfs(grafo.metadatos.nodos_iniciales[0].id, grafo.informacion_heuristica);
    } else {
    std::cerr << "[ERROR] grafo.metadatos.nodos_iniciales está vacío. Revisa el archivo de entrada o la inicialización." << std::endl;
    exit(1);  // o maneja el error apropiadamente
    }

    //bfs(grafo.metadatos.nodos_iniciales[0].id, grafo.informacion_heuristica);
    //bfs(32308, grafo.informacion_heuristica);
    //////////////////TEST

    //aco->abrir_file();
    auto start = std::chrono::high_resolution_clock::now();
    for (; aco->epoca_actual < aco->epocas; aco->epoca_actual++)
    {
        cout << "⌚"<< "Época " << aco->epoca_actual << endl;
        aco->resolver();
        aco->reset();
    }

    auto final = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_tiempo = final - comienzo;
    std::cout << "Tiempo de ejecución: " << total_tiempo.count() << " segundos" << std::endl;


    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //aco->cerrar_file();
    aco->mostrar_solucion(config.show_solucion);
    if (!parametros_base.silence || parametros_base.rescol)
    {
        aco->exportar_solucion(duration, parametros_base);
    }
    // if (!parametros_base.silence){
    //     aco->exportar_mapa_resultados();
    //     std::string archivo_salida = aco->get_filename();

    //     std::stringstream ss;
    //     ss << "python Grafico.py " << archivo_salida << " " << config.show_grafico;
    //     std::string comando = ss.str();

    //     std::stringstream ss2;
    //     ss2 << "python Visualizador.py " << parametros_base.nombre_instancia << " " << archivo_salida << " " << config.show_grafico;
    //     std::string comando2 = ss2.str();
        
    //     cout << "Programa finalizado correctamente" << endl;
    //     for (int i = 0; i < parametros_base.num_hormigas; i++)
    //         cout << "🐜 ";
    //     cout << endl;
    //     cout << endl;
        

    //     std::system(comando.c_str());    
    //     std::system(comando2.c_str());     
    // }

    //liberar_memoria(grafo);
    delete aco;
    return 0;
}
