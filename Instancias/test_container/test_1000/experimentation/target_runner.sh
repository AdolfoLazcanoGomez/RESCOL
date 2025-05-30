#!/bin/bash


CONFIG_ID="$1"
INSTANCE_ID="$2"
SEED="$3"
INSTANCE="$4"
shift 5
CONFIG_PARAMS=$*



# Definir variables para los parámetros restantes
metodo=0
tiempo_max=600
epocas=1
beta0="--beta0"
silence="--silence"
dir_salida="./resultados"
prefijo_salida="nuevos"
rescol="--rescol"
# Crear un identificador único para el archivo de registro
timestamp=$(date +%Y%m%d%H%M%S%N)
log_file="rescol_log_$timestamp.txt"



# Ejecutar el comando con los parámetros y redirigir la salida al archivo de registro
/home/alazcano/RESCOL/RESCOL $INSTANCE --semilla $SEED ${CONFIG_PARAMS} --salida-dijkstra --usar-tiempo --tiempo-max $tiempo_max --metodo $metodo  --epocas $epocas  $beta0  $usar_limitador $silence --dir-salida $dir_salida --prefijo-salida $prefijo_salida $rescol &> $log_file

# Extraer el número después de "Mejor costo:"
mejor_costo=$(grep -oP 'Mejor costo: \K[0-9.]*' $log_file)

echo "$mejor_costo"
rm $log_file