# Definir variables
#EXECUTABLE="/home/alazcano/Rescol_NM/Rescol_NM/RESCOL"
#INPUT_FILE="/home/alazcano/Rescol_NM/Rescol_NM/Instancias/CasoRealMediano.txt"
#INSTANCES_DIR="/home/alazcano/Rescol_NM/Rescol_NM/Instancias/Instancias_test"
#INSTANCE_NAME=$(basename "$INPUT_FILE" .txt) # Nombre de la instancia sin extensiÃ³n

# Obtener rutas relativas al directorio del script
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")
EXECUTABLE="$PROJECT_ROOT/RESCOL"
INSTANCES_DIR="$PROJECT_ROOT/Instancias/CN"
MULTIPLICITY_DIR="$SCRIPT_DIR/multiplicidad"

METODO="0"
#ITER_MAX="1000"
NUM_HORMIGAS="26"
EPOCAS="1"
salida_Floyd_Warshall="salida-Floyd-Warshall"

BETA0="--beta0"
ALFA="4.99"
RHO="0.3"
TAU="2.06"

#USAR_ITERACIONES="--usar-iteraciones"
#USAR_LIMITADOR="--usar-limitador"
USAR_TIEMPO="--usar-tiempo"
TIEMPO_MAX="30"
#VALOR_LIMITADOR="1"
SILENCE="--silence"
#DIR_SALIDA="./resultados"
PREFIJO_SALIDA="nuevos"
RESCOL="--rescol"

## Ejecutar el comando 10 veces con la primera semilla fija y las siguientes secuenciales

RESULTS_ROOT="$SCRIPT_DIR/results"
mkdir -p "$RESULTS_ROOT"

#SEMILLA_INICIAL=1234567

for INPUT_FILE in "$INSTANCES_DIR"/*.txt; do
  SEMILLA_INICIAL=1234567
  if [[ ! -f "$INPUT_FILE" ]]; then
    echo "No se encontraron archivos .txt en $INSTANCES_DIR"
    break
  fi

  INSTANCE_NAME=$(basename "$INPUT_FILE" .txt)
  SEMILLA=$SEMILLA_INICIAL

  echo "Procesando instancia: $INSTANCE_NAME"

  INSTANCE_RESULTS_DIR="$RESULTS_ROOT/$INSTANCE_NAME"
  mkdir -p "$INSTANCE_RESULTS_DIR"

  MULTIPLICITY_FILE="$MULTIPLICITY_DIR/${INSTANCE_NAME}.txt"
  LIMITADOR_ARGS=""
  if [[ -f "$MULTIPLICITY_FILE" ]]; then
    if VALOR_LIMITADOR=$(awk 'NR>2 && $3 ~ /^[0-9]+$/ { if ($3 > max) max = $3 } END { print (max ? max : 1) }' "$MULTIPLICITY_FILE"); then
      LIMITADOR_ARGS="--usar-limitador --valor-limitador $VALOR_LIMITADOR"
    else
      echo "Advertencia: no se pudo leer el archivo de multiplicidad para $INSTANCE_NAME" >&2
    fi
  else
    echo "Advertencia: no se encontró archivo de multiplicidad para $INSTANCE_NAME" >&2
  fi
  
  for i in {0..19}; do
    #OUTPUT_FILE="${INSTANCE_NAME}_$SEMILLA.txt"
    OUTPUT_FILE="$INSTANCE_RESULTS_DIR/${INSTANCE_NAME}_$SEMILLA.txt"
    DIR_SALIDA="$INSTANCE_RESULTS_DIR"
    echo "  Ejecutando con semilla: $SEMILLA, salida: $OUTPUT_FILE"
    # "$EXECUTABLE" "$INPUT_FILE" --metodo $METODO --alfa $ALFA --rho $RHO --tau-as $TAU \
    #   --tiempo-max $TIEMPO_MAX --num-hormigas $NUM_HORMIGAS --epocas $EPOCAS $salida_Floyd_Warshall \
    #   $BETA0 $USAR_TIEMPO $USAR_LIMITADOR --valor-limitador $VALOR_LIMITADOR $SILENCE \
    #   --semilla $SEMILLA --dir-salida "$DIR_SALIDA" --prefijo-salida $PREFIJO_SALIDA $RESCOL > "$OUTPUT_FILE"
    "$EXECUTABLE" "$INPUT_FILE" --metodo $METODO --alfa $ALFA --rho $RHO --tau-as $TAU \
      --tiempo-max $TIEMPO_MAX --num-hormigas $NUM_HORMIGAS --epocas $EPOCAS $salida_Floyd_Warshall \
      $BETA0 $USAR_TIEMPO $LIMITADOR_ARGS $SILENCE \
      --semilla $SEMILLA --dir-salida "$DIR_SALIDA" --prefijo-salida $PREFIJO_SALIDA $RESCOL > "$OUTPUT_FILE"
    SEMILLA=$((SEMILLA + 1))
  done
done