#ifndef OSCILADOR_H
#define OSCILADOR_H

#include <vector>


struct ParametroOscilante {
    float& valor;
    double incremento;
    double limiteInferior;
    double limiteSuperior;
    bool incrementando;

    ParametroOscilante(float& ref, double inc, double min, double max);

    void oscilar();
    void oscilar_caotico();
};


class Oscilador {
private:
    std::vector<ParametroOscilante> parametros;

public:
    Oscilador();

    void agregarParametro(float& ref, double inc, double min, double max);

    void oscilar();
    void oscilar_caotico();
};

#endif
