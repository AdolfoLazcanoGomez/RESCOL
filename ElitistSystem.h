#ifndef ELITS_H
#define ELITS_H

#include "aco.h"

class ElitistSystem : public ACO
{
protected:

    struct ArcoInfo {
        double calidad;
        bool visitado;
        int vecesVisitado;
    };
    ACOArgs parametros;    
    //double q_0; // umbral de probabilidad ACS
    double q;
    //double csi;
    std::unordered_map<int, std::vector<std::pair<Arco *, ArcoInfo>>> mapaACS; 


public:
    ElitistSystem(Graph *instancia, ACOArgs parametros);
    void resolver() override;
    void iterar() override;
    void inicializar_feromonas() override;
    void set_parametrosACS(ACOArgs parametros);
private:
};

#endif
