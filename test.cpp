#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include "graph.h"


void bfs(int nodoInicial, std::unordered_map<int, std::unordered_map<Arco *, double>> &adyacencia) {
    std::queue<int> cola;
    std::unordered_map<int, bool> visitado;

    // Marcar todos los nodos como no visitados al inicio
    for (auto &par : adyacencia) {
        visitado[par.first] = false;
    }

    // Marcar el nodo inicial como visitado y encolarlo
    visitado[nodoInicial] = true;
    cola.push(nodoInicial);

    while (!cola.empty()) {
        int nodoActual = cola.front();
        cola.pop();

        // Procesar los vecinos del nodo actual
        for (auto &par : adyacencia[nodoActual]) {
            Arco *arco = par.first;
            int idVecino = arco->destino->id;

            if (!visitado[idVecino]) {
                visitado[idVecino] = true;
                cola.push(idVecino);
            }
        }
    }
    /*
    int contador = 0;
    for (auto &par : adyacencia) {
        if (visitado[par.first] == false){
            std::cout << "El nodo " << par.first << " no fue visitado" << std::endl;
            contador++;
        }
    }
    cout << "cantidad de nodos no visitados: " << contador << endl;
    contador = 0;
    for (auto &par : adyacencia) {
        if (visitado[par.first] != false){
            std::cout << "El nodo " << par.first << " fue visitado" << std::endl;
            contador++;
        }
    }
    cout << "cantidad de nodos visitados: " << contador << endl;
    */
}


    
