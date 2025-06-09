#include "graph.h"

Graph::Graph(){
    for (auto &par : arcos){
        delete par.second;
    }
    arcps.clear();
}