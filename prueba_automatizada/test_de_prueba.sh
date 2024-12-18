#!/bin/bash

# Directorio de las instancias
instancias_dir="/home/pcontreras/experimentos_rescol/Instancias_training"

# Comando base
comando_base="/home/pcontreras/rescolants/RESCOL"

# Parámetros fijos
parametros="--metodo 0 --num-hormigas 30 --salida-dijkstra --beta0 --alfa 4.81 --rho 0.24 --tau-as 3.56 --valor-limitador 3 --epocas 1 --valor-sin-nuevas-visitas 1 --usar-tiempo --tiempo-max 1 --silence false --dir-salida ./resultados --prefijo-salida nuevos --rescol --semilla 1"

# Contador de instancias
contador_instancias=0
contador_exitos=0

# Iterar sobre cada archivo de instancia en el directorio
for instancia in "$instancias_dir"/*.txt; do
    # Incrementar contador de instancias
    ((contador_instancias++))
    
    # Construir el comando completo
    nombre_instancia=$(basename "$instancia" .txt)
    log_file="/home/pcontreras/rescolants/prueba_automatizada/${nombre_instancia}_log.txt"

    comando_completo="$comando_base $instancia $parametros"
    echo "Ejecutando $nombre_instancia"
    # Ejecutar el comando
    $comando_completo >> "$log_file"
    
    # Verificar si el comando se ejecutó correctamente
    if [ $? -ne 0 ]; then
        echo "Error ejecutando el comando para la instancia: $instancia"
        exit 1
    fi
    
    # Verificar la aparición de "Mejor costo: " en el archivo de log
    mejor_costo_count=$(grep -c "Mejor costo: " "$log_file")
    
    if [ "$mejor_costo_count" -gt 0 ]; then
        ((contador_exitos++))
    else
        echo "Error: No se encontró 'Mejor costo: ' en el log de $nombre_instancia"
    fi
done

# Verificar si la cantidad de éxitos es igual al número de instancias
if [ "$contador_exitos" -eq "$contador_instancias" ]; then
    echo "Prueba exitosa sin fallos"
else
    echo "Prueba con fallos"
fi

rm *_log.txt
