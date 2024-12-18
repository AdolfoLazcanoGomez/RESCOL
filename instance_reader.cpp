#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <set>

#include "instance_reader.h"
#include "graph.h"
#include "helpers.h"

using namespace std;

// Función para leer una instancia
/* Parámetros:
    - nombre_archivo: Nombre del archivo que contiene la instancia
    - leer_restricciones: Indica si se leerán las restricciones del archivo
    - leer_coordenadas: Indica si se leerán las coordenadas del archivo
*/
Graph leerInstancia(const std::string &nombre_archivo,
                    bool leer_restricciones,
                    bool leer_coordenadas)
{
    // Inicializar el grafo
    Graph g;
    std::unordered_map<std::string, std::string> encabezado;
    std::string lineaDato;
    std::string ruta_archivo;
    std::set<std::pair<int, int>> arcosCreados;
    std::set<std::pair<int, int>> arcosCreadosOpcionales;
    int IdArco = 0;
    int id_secundario = 1;

    ruta_archivo = nombre_archivo;
    std::ifstream infile(ruta_archivo);
    if (!infile)
    {
        std::cerr << "No se pudo abrir el archivo: " << nombre_archivo << std::endl;
        return Graph();
    }

    // Leer la información general del archivo
    for (int i = 0; i < 8; i++)
    {
        std::string linea;
        std::getline(infile, linea);
        size_t pos = linea.find(':');
        std::string clave = eliminarEspacios(linea.substr(0, pos));
        std::string valor = eliminarEspacios(linea.substr(pos + 1));
        encabezado[clave] = valor;
    }    
    g.metadatos.encabezado = encabezado;

    // Saltarse header aristas obligatorias
    std::getline(infile, lineaDato);

    // Leer aristas obligatorias
    for (int i = 0; i < std::stoi(encabezado["ARISTAS_REQ"]); i++)
    {
        id_secundario = 1;
        std::getline(infile, lineaDato);
        std::istringstream arcStream(lineaDato);
        std::string bi_or_uni;
        double id, origen, destino, costo_recorrido, costo_recoleccion;
        if (!(arcStream >> bi_or_uni >> origen >> destino >> costo_recorrido >> costo_recoleccion))
        {
            // Error de formato
            std::cerr << "Error en el formato del archivo" << std::endl;
            return Graph();
        }

        // Crear los nodos si no existen
        if (g.nodos.count(origen) == 0)
        {
            Nodo nodo_origen;
            nodo_origen.id = origen;
            g.nodos[origen] = nodo_origen;
        }
        if (g.nodos.count(destino) == 0)
        {
            Nodo nodo_destino;
            nodo_destino.id = destino;
            g.nodos[destino] = nodo_destino;
        }

        // Crear el arco
        if (bi_or_uni == "uni" )
        {
            Arco *arco = new Arco;
            arco->id = IdArco;
            arco->costo_recorrido = costo_recorrido;
            arco->costo_recoleccion = costo_recoleccion;
            arco->origen = &g.nodos[origen];
            arco->destino = &g.nodos[destino];
            arco->bidireccional = false;
            g.arcos[IdArco] = arco;

            // Conectar el arco a sus nodos
            g.nodos[origen].saliente.push_back(*arco);
            g.nodos[destino].entrante.push_back(*arco);

            // nodos fantasma nodo origen
            Subnodo *subnodo_saliente_origen = new Subnodo;

            subnodo_saliente_origen->tipo = 2; // 2: saliente

            g.nodos[origen].subnodos.push_back(*subnodo_saliente_origen);

            // nodos fantasma nodo destino
            Subnodo *subnodo_entrante_destino = new Subnodo;


            subnodo_entrante_destino->tipo = 1; // 1: entrante

            g.nodos[destino].subnodos.push_back(*subnodo_entrante_destino);

            // Actualizar la información heurística
            double costo = arco->costo_recorrido;
            g.informacion_heuristica[origen][arco] = 1.0 / costo;

            IdArco++;
        }
        else // arco bidireccional
        {
            if (arcosCreados.find(std::make_pair(origen, destino)) == arcosCreados.end() &&
                arcosCreados.find(std::make_pair(destino, origen)) == arcosCreados.end())
            {

                Arco *arcoIda = new Arco;
                Arco *arcoVuelta = new Arco;

                arcoIda->id = IdArco;
                g.arcos[IdArco] = arcoIda;
                IdArco++;

                arcoVuelta->id = IdArco;
                g.arcos[IdArco] = arcoVuelta;
                IdArco++;

                arcoIda->costo_recorrido = costo_recorrido;
                arcoIda->costo_recoleccion = costo_recoleccion;

                arcoVuelta->costo_recorrido = costo_recorrido;
                arcoVuelta->costo_recoleccion = costo_recoleccion;

                arcoIda->origen = &g.nodos[origen];
                arcoIda->destino = &g.nodos[destino];

                arcoVuelta->origen = &g.nodos[destino];
                arcoVuelta->destino = &g.nodos[origen];

                arcoIda->bidireccional = true;
                arcoVuelta->bidireccional = true;

                arcoIda->arco_reciproco = arcoVuelta;
                arcoVuelta->arco_reciproco = arcoIda;

                // Conectar el arco a sus nodos
                g.nodos[origen].saliente.push_back(*arcoIda);
                g.nodos[destino].entrante.push_back(*arcoIda);

                g.nodos[origen].entrante.push_back(*arcoVuelta);
                g.nodos[destino].saliente.push_back(*arcoVuelta);

                // nodos fantasma nodo origen
                Subnodo *subnodo_entrante_origen = new Subnodo;
                Subnodo *subnodo_saliente_origen = new Subnodo;

                subnodo_entrante_origen->id_secundario = id_secundario;
                subnodo_saliente_origen->id_secundario = id_secundario;
                id_secundario++;

                subnodo_entrante_origen->tipo = 1; // 1: entrante
                subnodo_saliente_origen->tipo = 2; // 2: saliente

                subnodo_entrante_origen->nodo_reciproco = subnodo_saliente_origen;
                subnodo_saliente_origen->nodo_reciproco = subnodo_entrante_origen;

                g.nodos[origen].subnodos.push_back(*subnodo_entrante_origen);
                g.nodos[origen].subnodos.push_back(*subnodo_saliente_origen);

                // nodos fantasma nodo destino
                Subnodo *subnodo_entrante_destino = new Subnodo;
                Subnodo *subnodo_saliente_destino = new Subnodo;

                subnodo_entrante_destino->id_secundario = id_secundario;
                subnodo_saliente_destino->id_secundario = id_secundario;
                id_secundario++;

                subnodo_entrante_destino->tipo = 1; // 1: entrante
                subnodo_saliente_destino->tipo = 2; // 2: saliente

                subnodo_entrante_destino->nodo_reciproco = subnodo_saliente_destino;
                subnodo_saliente_destino->nodo_reciproco = subnodo_entrante_destino;

                g.nodos[destino].subnodos.push_back(*subnodo_entrante_destino);
                g.nodos[destino].subnodos.push_back(*subnodo_saliente_destino);

                // Actualizar la información heurística
                double costo = arcoIda->costo_recorrido;
                g.informacion_heuristica[origen][arcoIda] = 1.0 / costo;

                costo = arcoVuelta->costo_recorrido;
                g.informacion_heuristica[destino][arcoVuelta] = 1.0 / costo;

                arcosCreados.insert(std::make_pair(origen, destino));
                arcosCreados.insert(std::make_pair(destino, origen));
            }
        }
    }

    // Saltarse header aristas opcionales
    std::getline(infile, lineaDato);

    // Leer aristas opcionales
    if (encabezado["ARISTAS_NOREQ"] != "0")
    {
     for (int i = 0; i < std::stoi(encabezado["ARISTAS_NOREQ"]); i++)
    {
        id_secundario = 1;
        std::getline(infile, lineaDato);
        std::istringstream arcStream(lineaDato);
        std::string bi_or_uni;
        double id, origen, destino, costo_recorrido, costo_recoleccion;
        if (!(arcStream >> bi_or_uni >> origen >> destino >> costo_recorrido >> costo_recoleccion))
        {
            // Error de formato
            std::cerr << "Error en el formato del archivo" << std::endl;
            return Graph();
        }

        // Crear los nodos si no existen
        if (g.nodos.count(origen) == 0)
        {
            Nodo nodo_origen;
            nodo_origen.id = origen;
            g.nodos[origen] = nodo_origen;
        }
        if (g.nodos.count(destino) == 0)
        {
            Nodo nodo_destino;
            nodo_destino.id = destino;
            g.nodos[destino] = nodo_destino;
        }

        // Crear el arco
        if (bi_or_uni == "uni" )
        {
            Arco *arco = new Arco;
            arco->id = IdArco;
            arco->obligatoria = false;
            arco->costo_recorrido = costo_recorrido;
            arco->costo_recoleccion = costo_recoleccion;
            arco->origen = &g.nodos[origen];
            arco->destino = &g.nodos[destino];
            arco->bidireccional = false;
            g.arcos[IdArco] = arco;

            // Conectar el arco a sus nodos
            g.nodos[origen].saliente.push_back(*arco);
            g.nodos[destino].entrante.push_back(*arco);

            // nodos fantasma nodo origen
            Subnodo *subnodo_saliente_origen = new Subnodo;

            subnodo_saliente_origen->tipo = 2; // 2: saliente

            g.nodos[origen].subnodos.push_back(*subnodo_saliente_origen);

            // nodos fantasma nodo destino
            Subnodo *subnodo_entrante_destino = new Subnodo;


            subnodo_entrante_destino->tipo = 1; // 1: entrante

            g.nodos[destino].subnodos.push_back(*subnodo_entrante_destino);

            // Actualizar la información heurística
            double costo = arco->costo_recorrido;
            g.informacion_heuristica[origen][arco] = 1.0 / costo;

            IdArco++;
        }
        else // arco bidireccional
        {
            if (arcosCreados.find(std::make_pair(origen, destino)) == arcosCreados.end() &&
                arcosCreados.find(std::make_pair(destino, origen)) == arcosCreados.end())
            {

                Arco *arcoIda = new Arco;
                Arco *arcoVuelta = new Arco;

                arcoIda->id = IdArco;
                arcoIda->obligatoria = false;
                g.arcos[IdArco] = arcoIda;
                IdArco++;

                arcoVuelta->id = IdArco;
                g.arcos[IdArco] = arcoVuelta;
                IdArco++;

                arcoIda->costo_recorrido = costo_recorrido;
                arcoIda->costo_recoleccion = costo_recoleccion;

                arcoVuelta->costo_recorrido = costo_recorrido;
                arcoVuelta->costo_recoleccion = costo_recoleccion;

                arcoIda->origen = &g.nodos[origen];
                arcoIda->destino = &g.nodos[destino];

                arcoVuelta->origen = &g.nodos[destino];
                arcoVuelta->destino = &g.nodos[origen];

                arcoIda->bidireccional = true;
                arcoVuelta->bidireccional = true;

                arcoIda->arco_reciproco = arcoVuelta;
                arcoVuelta->arco_reciproco = arcoIda;

                // Conectar el arco a sus nodos
                g.nodos[origen].saliente.push_back(*arcoIda);
                g.nodos[destino].entrante.push_back(*arcoIda);

                g.nodos[origen].entrante.push_back(*arcoVuelta);
                g.nodos[destino].saliente.push_back(*arcoVuelta);

                // nodos fantasma nodo origen
                Subnodo *subnodo_entrante_origen = new Subnodo;
                Subnodo *subnodo_saliente_origen = new Subnodo;

                subnodo_entrante_origen->id_secundario = id_secundario;
                subnodo_saliente_origen->id_secundario = id_secundario;
                id_secundario++;

                subnodo_entrante_origen->tipo = 1; // 1: entrante
                subnodo_saliente_origen->tipo = 2; // 2: saliente

                subnodo_entrante_origen->nodo_reciproco = subnodo_saliente_origen;
                subnodo_saliente_origen->nodo_reciproco = subnodo_entrante_origen;

                g.nodos[origen].subnodos.push_back(*subnodo_entrante_origen);
                g.nodos[origen].subnodos.push_back(*subnodo_saliente_origen);

                // nodos fantasma nodo destino
                Subnodo *subnodo_entrante_destino = new Subnodo;
                Subnodo *subnodo_saliente_destino = new Subnodo;

                subnodo_entrante_destino->id_secundario = id_secundario;
                subnodo_saliente_destino->id_secundario = id_secundario;
                id_secundario++;

                subnodo_entrante_destino->tipo = 1; // 1: entrante
                subnodo_saliente_destino->tipo = 2; // 2: saliente

                subnodo_entrante_destino->nodo_reciproco = subnodo_saliente_destino;
                subnodo_saliente_destino->nodo_reciproco = subnodo_entrante_destino;

                g.nodos[destino].subnodos.push_back(*subnodo_entrante_destino);
                g.nodos[destino].subnodos.push_back(*subnodo_saliente_destino);

                // Actualizar la información heurística
                double costo = arcoIda->costo_recorrido;
                g.informacion_heuristica[origen][arcoIda] = 1.0 / costo;

                costo = arcoVuelta->costo_recorrido;
                g.informacion_heuristica[destino][arcoVuelta] = 1.0 / costo;

                arcosCreados.insert(std::make_pair(origen, destino));
                arcosCreados.insert(std::make_pair(destino, origen));
            }
        }
    }
    }

    // Saltarse header Coordenadas
    std::getline(infile, lineaDato);

    // Leer coordenadas
    // TODO: Formato coordenadas
    // Datos dummy ⚠️
    if (leer_coordenadas)
    {
        for (int i = 0; i < std::stoi(encabezado["VERTICES"]); i++)
        {
            std::getline(infile, lineaDato);
            std::istringstream arcStream(lineaDato);
            int id;
            double x, y;
            if (!(arcStream >> x >> y))
            {
                // Error de formato
                std::cerr << "Error en el formato del archivo" << std::endl;
                return Graph();
            }
        }
    }
    else
    {
        /*
        Se elimina en siguientes versiones
        for (int i = 0; i < std::stoi(encabezado["VERTICES"]); i++)
            std::getline(infile, lineaDato);
        */          
    }

    // Saltarse header Restricciones
    std::getline(infile, lineaDato);

    // Leer restricciones
    // TODO: Evaluar si son necesarias las restricciones
    if (encabezado["RESTRICCIONES"] != "0")
    {
        if (leer_restricciones)
        {
            for (int i = 0; i < std::stoi(encabezado["RESTRICCIONES"]); i++)
            {
                std::getline(infile, lineaDato);
                std::istringstream arcStream(lineaDato);
                int id, nodo1, nodo2, nodo3;
                if (!(arcStream >> nodo1 >> nodo2 >> nodo3))
                {
                    // Error de formato
                    std::cerr << "Error en el formato del archivo" << std::endl;
                    return Graph();
                }
            }
        }
        else
        {
            for (int i = 0; i < std::stoi(encabezado["RESTRICCIONES"]); i++)
            {
                std::getline(infile, lineaDato);
            }
        }
    }

    // Saltarse header Depositos o nodos iniciales
    std::getline(infile, lineaDato);

    // Leer depositos o nodos iniciales
    for (int i = 0; i < std::stoi(encabezado["NODOS_INICIALES"]); i++)
    {
        std::getline(infile, lineaDato);
        std::istringstream arcStream(lineaDato);
        int id, nodo;
        if (!(arcStream >> nodo))
        {
            // Error de formato
            std::cerr << "Error en el formato del archivo" << std::endl;
            return Graph();
        }
        g.metadatos.nodos_iniciales.push_back(g.nodos[nodo]);
    }

    // Saltarse header Terminal
    std::getline(infile, lineaDato);

    // Leer terminales
    if (encabezado["NODOS_TERMINO"] != "0")
    {
        for (int i = 0; i < std::stoi(encabezado["NODOS_TERMINO"]); i++)
        {
            std::getline(infile, lineaDato);
            std::istringstream arcStream(lineaDato);
            int id, nodo;
            if (!(arcStream >> nodo))
            {
                // Error de formato
                std::cerr << "Error en el formato del archivo" << std::endl;
                return Graph();
            }
            g.metadatos.nodos_termino.push_back(g.nodos[nodo]);
        }
    }
    else
    {
        for (int i = 0; i < std::stoi(encabezado["VERTICES"]); i++)
        {
            g.metadatos.nodos_termino.push_back(g.nodos[i + 1]);
        }
    }

    // Cerrar el archivo
    infile.close();

    // Crear subnodos por cada nodo

    for (auto &par : g.nodos)
    {
        for (auto &par2 : par.second.subnodos)
        {
            if (par2.tipo == 1) // 1:entrante
            {
                for (auto &par3 : par.second.subnodos)
                {
                    if ((par3.tipo == 2) && (par3.id_secundario != par2.id_secundario)) // 2:saliente
                    {
                        Arco *arco = new Arco;
                        arco->id = IdArco;
                        IdArco++;

                        par2.saliente.push_back(*arco);
                        par3.entrante.push_back(*arco);
                    }
                }
            }
        }
        
        
    }

    return g;
}