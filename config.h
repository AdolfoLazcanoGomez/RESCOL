#ifndef CONFIG_H
#define CONFIG_H

#include "argparse/argparse.hpp"
#include "enums.h"

struct ACOArgs : public argparse::Args
{
    // IN/OUT
    std::string &nombre_instancia = arg("Nombre de la instancia").set_default("CasoRealChiquito.txt");
    //bool &crear_directorios_solucion = flag("crear-directorios-solucion", "permite crear directorios para las soluciones");
    std::string &directorio_salida = kwarg("dir-salida", "Directorio para archivos de salida").set_default("./resultados");
    std::string &prefijo_salida = kwarg("prefijo-salida", "Prefijo para archivos de salida").set_default("resultados");
    bool &rescol = flag("rescol", "formato entrada y salida para la interfaz rescol");
    bool &irace = flag("irace", "usa el formato de llamada e instancias de irace");
    bool &silence = flag("silence", "output minimo");
    bool &conectar_bd = flag("bdconn", "Flag que permite o no subir los resultados a la BD");
    
    // Semilla
    double &semilla = kwarg("semilla", "semilla de generacion aleatoria").set_default(1188444438);
    bool &full_aleatorio = flag("random-seed", "semilla aleatoria");
    
    // Budget
    int &iteraciones_max = kwarg("iter-max", "Cantidad de iteraciones máximas").set_default(1000);
    int &evaluaciones_maximas = kwarg("eval-max", "Cantidad de evaluaciones máximas").set_default(500);
    double &tiempo_maximo = kwarg("tiempo-max", "Cantidad de tiempo máximo").set_default(30.0);
    bool &usar_iteraciones = flag("usar-iteraciones", "usa iteraciones en vez de evaluaciones o tiempo,  no se pueden usar mas de una a la vez");
    bool &usar_evaluaciones = flag("usar-evaluaciones", "usa evaluaciones en vez de iteraciones o tiempo,  no se pueden usar mas de una a la vez");
    bool &usar_tiempo = flag("usar-tiempo", "usa tiempo en vez de iteraciones o evaluaciones, no se pueden usar mas de una a la vez");
    int &epocas = kwarg("epocas", "Épocas").set_default(3);
    
    // AS basico
    int &metodo = kwarg("metodo", "Metodo de resolucion").set_default(1);
    int &num_hormigas = kwarg("num-hormigas", "Número de hormigas").set_default(10);
    double &umbral_inferior = kwarg("umbral-inf", "Umbral inferior para las feromonas").set_default(0.001);
    float &tau = kwarg("tau-as", "Parámetro tau, asociado a las feromonas iniciales").set_default(1.0);
    float &alfa = kwarg("alfa", "Parámetro alfa").set_default(1.0);
    float &beta = kwarg("beta", "Parámetro beta").set_default(2.0);
    bool &beta_0 = flag("beta0", "Flag que fuerza el valor de beta a 0 en el recorrido");
    float &rho = kwarg("rho", "Parámetro rho, asociado a la evaporación de feromonas").set_default(0.5);
    
    // Matriz de salida
    float &rho_secundario = kwarg("rho-sec", "Parámetro rho, asociado a la evaporación de feromonas").set_default(0.8);
    bool &usaMatrizSecundaria = flag("matriz-secundaria", "Flag que permite o no el uso de la matriz de salida");
    float &beta_salida = kwarg("beta-salida", "Parámetro beta matriz de salida").set_default(2.0);
    float &rho_salida = kwarg("rho-salida", "Parámetro rho matriz de salida, asociado a la evaporación de feromonas").set_default(0.8);

    // Algoritmo salida
    bool &usaDijkstra = flag("salida-dijkstra", "Flag que permite o no el uso del algoritmo de Dijkstra");

    // Limitador
    bool &limitador = flag("usar-limitador", "usa el limitador de pasadas");
    int &valor_limitador = kwarg("valor-limitador", "cantidad limite de veces que se puede pasar por un arco").set_default(4);
    
    // Uso de modo sin mejora de las calles
    //int &UsaLimiteDeMejoras = kwarg("usar-sin-nuevos-visitas", "Usar el limite sin nuevas calles").set_default(4);
    int &LimiteDeMejoras = kwarg("valor-sin-nuevas-visitas", "Cantidad de visitas sin nuevas calles").set_default(1);

    // MinMax
    double &umbral_superior = kwarg("umbral-sup", "Umbral superior para las feromonas").set_default(161161.161);
    int &umbral_sin_mejora_limite = kwarg("umbral-sin-mejora", "Cantidad de iteraciones sin mejora para la actualización de feromonas").set_default(10);
    int &a = kwarg("a-mm", "Parámetro a, asociado a la actualización de feromonas").set_default(13);


    // Elitist
    double &factor_elitista = kwarg("factor-elitist", "Parametro asociado a la estrategia elitista de feromonas").set_default(1.0);

    // ACS
    double &q_0 = kwarg("q0", "Parámetro q0").set_default(0.389); // umbral de probabilidad ACS
    double &csi = kwarg("csi", "Parámetro de actualizacion local feromonas").set_default(0.1);

    // Oscilador
    int &oscilador = kwarg("oscilador", "Parametro que permite o no el uso del oscilador de parametros, 0 es desactivado, 1 es stepper, 2 es caotico").set_default(0);
    float &inc_alfa = kwarg("inc-alfa", "Incremento de alfa").set_default(0.1);
    float &min_alfa = kwarg("min-alfa", "Valor minimo de alfa").set_default(1);
    float &max_alfa = kwarg("max-alfa", "Valor maximo de alfa").set_default(5);
    float &inc_beta = kwarg("inc-beta", "Incremento de beta").set_default(0.1);
    float &min_beta = kwarg("min-beta", "Valor minimo de beta").set_default(1);
    float &max_beta = kwarg("max-beta", "Valor maximo de beta").set_default(0.5);
    float &inc_rho = kwarg("inc-rho", "Incremento de rho").set_default(0.1);
    float &min_rho = kwarg("min-rho", "Valor minimo de rho").set_default(0.1);
    float &max_rho = kwarg("max-rho", "Valor maximo de rho").set_default(0.5);
    float &inc_rho_secundario = kwarg("inc-rho-sec", "Incremento de rho secundario").set_default(0.1);
    float &min_rho_secundario = kwarg("min-rho-sec", "Valor minimo de rho secundario").set_default(0.1);
    float &max_rho_secundario = kwarg("max-rho-sec", "Valor maximo de rho secundario").set_default(0.5);
    //float &inc_rho_salida = kwarg("inc-rho-salida", "Incremento de rho salida").set_default(0.1); No usar oscilador en mecanismos de salida
    //float &min_rho_salida = kwarg("min-rho-salida", "Valor minimo de rho salida").set_default(0.1);
    //float &max_rho_salida = kwarg("max-rho-salida", "Valor maximo de rho salida").set_default(0.5);


    
};


struct ConfigPrograma
{
    bool leer_restricciones = false;
    bool leer_coordenadas = false;
    bool debug = true;
    bool debug_ACO = true;
    bool show_solucion = true;    
    std::string show_grafico = "False";
    
};

#endif
