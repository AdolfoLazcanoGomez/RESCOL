#ifndef INSTANCE_READER_H
#define INSTANCE_READER_H

#include "graph.h" 
#include <string>
#include <vector>
using namespace std;
Graph leerInstancia(const std::string &filename,
                    bool leer_restricciones,
                    bool leer_coordenadas);

std::vector<std::vector<double>> floydWarshall(Graph &g);

#endif
