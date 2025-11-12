#ifndef ACO_H
#define ACO_H


#include <math.h>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <filesystem>
#include <chrono>

#include "graph.h"
#include "config.h"
#include "oscilador.h"

// Representación de las Feromonas
struct Feromona
{
    Nodo *nodo_inicial = nullptr;  // Nodo inicial
    Nodo *nodo_terminal = nullptr; // Nodo terminal
    double cantidad = 0.0;         // Cantidad de feromonas
    
};

struct Camion
{
    //NUEVO
    std::vector<Arco> camino_tour;                          // Camino recorrido hasta completar los requisitos (capacidad + tiempo)
    std::vector<Arco> camino_salida;                        // Camino de salida una vez completada la iteración
    std::vector<Arco> camino_final;                         // Suma de ambos caminos.
    double longitud_camino_tour = 0.0;
    double longitud_camino_salida = 0.0;
    double longitud_camino_final = 0.0;

    double capacidad_maxima = 900.0;
    double capacidad_restante = 900.0;
    double tiempo_disponible = 500.0;
    double tiempo_restante = 500.0;
    double costo_total = 0.0;
    double costo_camino_camion  = 0.0;
    double costo_total_camino_camion = 0.0;

};

// Representación de una Hormiga
struct Hormiga
{
    Nodo *nodo_actual = nullptr;                            // Nodo actual
    std::unordered_map<Arco *, int> arcos_visitados_tour;   // Aristas visitadas
    std::unordered_map<Arco *, int> arcos_visitados_salida; // Aristas visitadas
    std::unordered_map<Arco *, int> arcos_visitados_final;  // Aristas visitadas

    int saltosTour = 0;                                     // Cantidad de saltos que toma terminar un tour completamente
    int saltosSalida = 0;                                   // Cantidad de saltos desde el termino del tour hasta la salida
                                                            // Longitud del camino, cuantas aristas se ha recorrido
    double costo_camino = 0.0;                              // Costo del camino, costo asociado a la recolección y recorrido
    bool solucion_valida = true;                            
    double longitud_final_camiones = 0.0;
    int id = 0;                                             // Identificador de la hormiga
    std::unordered_map<Arco *, Feromona> feromonas_locales; // Feromonas locales de la hormiga
    std::unordered_map<Arco*, int> arcos_no_visitados;      // Aristas no visitadas

    //NUEVO
    std::vector<Camion *> vector_camiones;                  // Camiones de cada hormiga
    std::vector<Arco> copia_camino_tour;                  // Copia arcos recorridos por hormiga
    
    //No sé si utilizarlos también ya que cada camion tiene los suyos y las hormigas no lo utilizaran luego 
    // std::vector<Arco> copia_camino_salida;                // Copia arcos recorridos para salir
    // std::vector<Arco> copia_camino_final;                 // Copia suma final arcos recorridos

};

class ACO
{
protected:
    Oscilador oscilador;

    virtual void iterar();                                 // Itera el algoritmo
    virtual void inicializar_feromonas() = 0;              // Itera el algoritmo
    Hormiga guardar_mejor_solucion_iteracion();            // Guarda la mejor solución
    bool solucionCompleta(Hormiga &hormiga);               // Verifica si la solución es completa
    bool enNodoTerminal(Hormiga &hormiga);                 // Verifica si la hormiga esta en un nodo final
    Hormiga mejor_solucion;                                // Mejor solución
    void limpiar();                                        // Limpia la memoria y datos del algoritmo
    void limpiar_rastro();
    void limpiar_camiones( Hormiga &hormiga);                                 // Limpia la memoria y datos del algoritmo
    std::unordered_map<Arco *, Feromona> feromonas;        // Feromonas
    std::unordered_map<Arco *, Feromona> feromonas_salida; // Feromonas
    std::vector<Hormiga> hormigas;     
    std::string nombre_metodo;
    std::string nombre_instancia_salida;
    void construirSolucion(Hormiga &hormiga);              // Construye la solución para una hormiga    
    int evaluaciones = 0;                                  // Evaluaciones de la funcion objetivo
    virtual Nodo *eligeSiguiente(Hormiga &hormiga, int &aux);        // Elige el siguiente nodo
    virtual void buscarSalidaMatriz(Hormiga &hormiga, int &aux);           // busca la salida
    virtual void buscarDijkstra(Hormiga &hormiga, int &aux);
    void calcular_longitud_camiones (Hormiga &hormiga);    // Calcula la longitud total de la lista de camiones por hormiga
    void calcular_costo_camino_camion (Hormiga &hormiga);  // Calcula el costo total de la lista de camiones de la hormiga.
    double total_camino_salida_camion(Hormiga &hormiga);   // Calcula longitud total de camino salida de camiones en hormiga.
    double total_camino_tour_camion(Hormiga &hormiga);     // Calcula longitud total camino tour de los camiones por homriga.
    virtual void visitar(Hormiga &hormiga, Nodo *nodo, int &aux);    // Visita el nodo siguiente
    Nodo* nodoInicialAleatorio (Graph *instancia);  //obtner nodo inicial aleatorio
    int mejor_costo = std::numeric_limits<int>::max();     // Mejor costo de la iteracion actual
    Graph *grafo = nullptr;                                // Grafo
    char filename[100];
    std::ofstream file;
    std::string prefijo_salida;
    std::filesystem::path directorio_salida;
    void set_parametros(ACOArgs parametros_base);
    float timeout = 600000 ;
    
    bool timeout_flag = false;
    bool timeout_flag_global = false;

    std::chrono::duration<double>  tiempo_total_ejecucion; // Tiempo total que se ejecuta el algoritmo
    std::chrono::high_resolution_clock::time_point comienzo_ejecucion; // Variable que toma el tiempo en que comenzo la ejecucion
    bool usar_tiempo = false;

    int usar_oscilador = 0;
    int valor_limitador = 999999;
    bool usar_bd = false;
    bool usar_iteraciones = false;
    bool usar_evaluaciones = false;

    int cantidad_sin_nuevo_arco = 0;

public:
    virtual ~ACO();
    //virtual void ejecutar() = 0;
    std::string nombre_instancia; // Nombre de la instancia
    std::string ruta_archivo_salida_csv;
    std::string ruta_archivo_config_salida_csv;
    float alfa;                   // Parámetro alfa
    int metodo;                   // Método (0: antsystem, 1: minmax, 3: ACS)
    float beta;                   // Parámetro beta
    float beta_salida;                   // Parámetro beta
    float rho;                    // Parámetro rho, asociado a la evaporacion de feromonas
    float tau;                    // Parámetro rho, asociado a la evaporacion de feromonas
    float rho_secundario;         // Parámetro rho, asociado a la evaporacion de feromonas
    float rho_salida;             // Parámetro rho de salida, asociado a la evaporacion de feromonas de la matriz de salida, deberia se un valor muy pequeño
    int iteraciones = 0;          // Cantidad de iteraciones, asociada a la funcion resolver y al criterio de parada
    int iteraciones_max;          // Cantidad de iteraciones maximas ||  = numevaluaciones / hormigas;    
    int evaluaciones_maximas;    // Cantidad de evaluaciones maximas
    bool debug = false;        // Flag que muestra o no informacion de debug como los caminos de las hormigas
    double umbral_inferior;    // Umbral inferior para las feromonas
    int num_hormigas;          // Numero de hormigas
    int epocas;                // Numero de epocas
    int epoca_actual = 0;      // Numero de epocas
    int LimiteDeMejoras;       // Cantidad limite que no se han elegido nuevos enlaces
    bool usarMatrizSecundaria; // Flag que controla el uso general de la matriz de salida, se pasa por parametros
    bool usaDijkstra;
    int sin_nuevo = 0;
    float acumulador_tiempo = 0;

    std::vector<std::vector<double>> dist;
    std::vector<std::vector<int>> next;
    

    ACO(Graph *instancia, ACOArgs parametros_base); // Constructor;
    // ACO(Graph *instancia, ParametrosACOBase parametros_base); // Constructor;
    void reset();

    virtual void resolver() = 0; // Resuelve el problema
    int contarArcosNoVisitadosObligatorios(Hormiga &hormiga);
    void mostrar_solucion(bool show_solucion); // Muestra la solución
    void set_mejor_feromonas();                // Setea las feromonas segun #hormigas/longitud_mejor_camino
    //void abrir_file();
    //void cerrar_file();
    std::string get_filename();
    //void set_filename(std::string filename);
    Hormiga get_mejor_solucion();
    void exportar_solucion(std::chrono::microseconds duration, ACOArgs parametros_base);
    void exportar_mapa_resultados();

    std::vector<Nodo*> findPath(Hormiga &hormiga, Nodo* targetNode);
    //std::vector<Nodo*> buscarRutaHastaArco(Hormiga &hormiga, Arco* arco, Graph &grafo);
    Nodo* ArcosNoVisitadoObligatorios(Hormiga &hormiga, int flag);

private:
    Hormiga mejor_hormiga; // Mejor solución de la iteracion actual
    int mejor_longitud = std::numeric_limits<int>::max(); // Mejor costo de la iteracion actual
};

#endif
