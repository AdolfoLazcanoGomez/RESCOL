#ifndef ELITS2_H
#define ELITS2_H

#include "aco.h"


class ElitistSystem2 : public ACO
{
protected:

public:
    double elitist;
    ElitistSystem2(Graph *instancia, ACOArgs parametros);
    void resolver() override;
    void iterar() override;
    void inicializar_feromonas() override;
    void set_parametrosACS(ACOArgs parametros);

};

#endif