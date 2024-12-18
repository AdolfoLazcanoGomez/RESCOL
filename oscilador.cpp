#include "oscilador.h"
#include "helpers.h"

ParametroOscilante::ParametroOscilante(float& ref, double inc, double min, double max)
    : valor(ref), incremento(inc), limiteInferior(min),
      limiteSuperior(max), incrementando(true) {}

void ParametroOscilante::oscilar() {
    if (incrementando) {
        valor += incremento;
        if (valor >= limiteSuperior) {
            valor = limiteSuperior;
            incrementando = false;
        }
    } else {
        valor -= incremento;
        if (valor <= limiteInferior) {
            valor = limiteInferior;
            incrementando = true;
        }
    }
}

void ParametroOscilante::oscilar_caotico() {
    valor = generar_numero_aleatorio(limiteInferior,limiteSuperior);
}

Oscilador::Oscilador() {}

void Oscilador::agregarParametro(float& ref, double inc, double min, double max) {
    parametros.emplace_back(ref, inc, min, max);
}

void Oscilador::oscilar() {
    for (auto& param : parametros) {
        param.oscilar();
    }
}
void Oscilador::oscilar_caotico() {
    for (auto& param : parametros) {
        param.oscilar_caotico();
    }
}
