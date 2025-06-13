#!/bin/bash

CONFIG_ID="$1"
INSTANCE_ID="$2"
SEED="$3"
INSTANCE="$4"
shift 5
CONFIG_PARAMS=$*

# Parámetros fijos del algoritmo
metodo=0
tiempo_max=480
epocas=1
beta0="--beta0"
silence="--silence"
dir_salida="./resultados"
prefijo_salida="nuevos"
rescol="--rescol"

# Identificadores únicos para los archivos temporales
timestamp=$(date +%Y%m%d%H%M%S%N)
log_file="rescol_log_${CONFIG_ID}_${timestamp}.txt"
err_file="rescol_err_${CONFIG_ID}_${timestamp}.txt"
cmd_file="rescol_cmd_${CONFIG_ID}_${timestamp}.txt"

# Registrar comando ejecutado
echo "=== EJECUCIÓN ===" >> "$log_file"
echo "Fecha y hora: $(date)" >> "$log_file"
echo "Semilla: $SEED" >> "$log_file"
echo "Instancia: $INSTANCE" >> "$log_file"
echo "Config ID: $CONFIG_ID" >> "$log_file"
echo "Instance ID: $INSTANCE_ID" >> "$log_file"
echo "Config Params: $CONFIG_PARAMS" >> "$log_file"
echo "==================" >> "$log_file"

# Registrar el comando exacto para replicación
echo "/home/alazcano/RESCOL/RESCOL \"$INSTANCE\" --semilla \"$SEED\" $CONFIG_PARAMS --salida-dijkstra --usar-tiempo --tiempo-max $tiempo_max --metodo $metodo --epocas $epocas $beta0 $usar_limitador $silence --dir-salida $dir_salida --prefijo-salida $prefijo_salida $rescol" > "$cmd_file"

# Ejecutar el algoritmo y capturar salida + estadísticas de memoria
/usr/bin/time -v /home/alazcano/RESCOL/RESCOL "$INSTANCE" --semilla "$SEED" $CONFIG_PARAMS \
  --salida-dijkstra --usar-tiempo --tiempo-max $tiempo_max --metodo $metodo --epocas $epocas \
  $beta0 $usar_limitador $silence --dir-salida $dir_salida --prefijo-salida $prefijo_salida $rescol \
  > "$log_file" 2> "$err_file"

exit_code=$?

# Detectar errores comunes
if [ $exit_code -eq 134 ]; then
  echo "999999 # ERROR: Abortado (posible double free, invalid free, etc.)"
elif [ $exit_code -eq 137 ]; then
  echo "999999 # ERROR: Finalizado por OOM (memoria insuficiente)"
elif [ $exit_code -ne 0 ]; then
  echo "999999 # ERROR: Código de salida $exit_code"
else
  # Extraer el mejor costo
  mejor_costo=$(grep -oP 'Mejor costo: \K[0-9.]+' "$log_file")
  max_mem=$(grep "Maximum resident set size" "$err_file" | awk '{print $6}')

  if [[ -z "$mejor_costo" ]]; then
    echo "999999 # ERROR: No se encontró 'Mejor costo' en el log - Memoria usada: ${max_mem:-N/A} KB"
  else
    echo "$mejor_costo # Memoria usada: ${max_mem:-N/A} KB"
  fi
fis

# Limpieza (deja cmd_file solo si hubo error)
if [ $exit_code -eq 0 ] && [[ -n "$mejor_costo" ]]; then
  rm -f "$log_file" "$err_file" "$cmd_file"
fi

#nohup /home/alazcano/R/x86_64-redhat-linux-gnu-library/4.0/irace/bin/irace --scenario scenario.txt --parallel 40 > irace.log 2>&1 &