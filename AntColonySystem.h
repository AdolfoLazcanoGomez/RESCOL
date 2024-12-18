#ifndef ACS_H
#define ACS_H

#include "aco.h"

class AntColonySystem : public ACO
{
protected:

    struct ArcoInfo {
        double calidad;
        bool visitado;
        int vecesVisitado;
    };
    ACOArgs parametros;    
    double q_0; // umbral de probabilidad ACS
    double q;
    double csi;
    std::unordered_map<int, std::vector<std::pair<Arco *, ArcoInfo>>> mapaACS; 


public:
    AntColonySystem(Graph *instancia, ACOArgs parametros);
    void resolver() override;
    void iterar() override;
    Nodo *eligeSiguiente(Hormiga &hormiga, int &aux); //override;
    //void visitar(Hormiga &hormiga, Nodo *nodo) override;
    void inicializar_feromonas() override;

    void set_parametrosACS(ACOArgs parametros);
private:
};

#endif
