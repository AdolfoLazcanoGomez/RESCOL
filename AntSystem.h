#ifndef AS_H
#define AS_H

#include "aco.h"


class AntSystem : public ACO
{
protected:

public:
    AntSystem(Graph *instancia, ACOArgs parametros);
    void resolver() override;
    void iterar() override;
    void inicializar_feromonas() override;
    void reemplazar_hormiga();

};

#endif
