# Definir variables
#EXECUTABLE="/home/alazcano/Rescol_NM/Rescol_NM/RESCOL"
#INPUT_FILE="/home/alazcano/Rescol_NM/Rescol_NM/Instancias/CasoRealMediano.txt"
#INSTANCES_DIR="/home/alazcano/Rescol_NM/Rescol_NM/Instancias/Instancias_test"
#INSTANCE_NAME=$(basename "$INPUT_FILE" .txt) # Nombre de la instancia sin extensiÃ³n

# Obtener rutas relativas al directorio del script
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
resolve_project_root() {
  local dir="$1"
  while [[ "$dir" != "/" ]]; do
    if [[ -d "$dir/.git" ]] || [[ -x "$dir/RESCOL" ]]; then
      printf '%s\n' "$dir"
      return 0
    fi
    dir=$(dirname "$dir")
  done
  return 1
}

PROJECT_ROOT=$(resolve_project_root "$SCRIPT_DIR")
if [[ -z "$PROJECT_ROOT" ]]; then
  echo "No se pudo localizar la raíz del proyecto desde $SCRIPT_DIR" >&2
  exit 1
fi

EXECUTABLE="$PROJECT_ROOT/RESCOL"
DEFAULT_INSTANCES_DIR="$PROJECT_ROOT/Instancias/QN"
INSTANCES_DIR="${INSTANCES_DIR:-$DEFAULT_INSTANCES_DIR}"
if [[ ! -d "$INSTANCES_DIR" ]]; then
  echo "No se encontró el directorio de instancias: $INSTANCES_DIR" >&2
  exit 1
fi

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
USAR_LIMITADOR="--usar-limitador"
USAR_TIEMPO="--usar-tiempo"
TIEMPO_MAX="30"
VALOR_LIMITADOR="1"
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
  INSTANCE_LOG="$INSTANCE_RESULTS_DIR/${INSTANCE_NAME}.log"
  : > "$INSTANCE_LOG"

  for i in {0..19}; do
    CURRENT_SEMILLA=$SEMILLA
    #OUTPUT_FILE="${INSTANCE_NAME}_$SEMILLA.txt"
    #OUTPUT_FILE="$INSTANCE_RESULTS_DIR/${INSTANCE_NAME}_$SEMILLA.txt"
    OUTPUT_FILE="$INSTANCE_RESULTS_DIR/${INSTANCE_NAME}_$CURRENT_SEMILLA.txt"
    DIR_SALIDA="$INSTANCE_RESULTS_DIR"
    echo "  Ejecutando con semilla: $CURRENT_SEMILLA, salida: $OUTPUT_FILE"
    {
      echo "=== Semilla $CURRENT_SEMILLA ==="
      "$EXECUTABLE" "$INPUT_FILE" --metodo $METODO --alfa $ALFA --rho $RHO --tau-as $TAU \
        --tiempo-max $TIEMPO_MAX --num-hormigas $NUM_HORMIGAS --epocas $EPOCAS $salida_Floyd_Warshall \
        $BETA0 $USAR_TIEMPO $USAR_LIMITADOR --salida-dijkstra --valor-limitador $VALOR_LIMITADOR $SILENCE \
        --semilla $CURRENT_SEMILLA --dir-salida "$DIR_SALIDA" --prefijo-salida $PREFIJO_SALIDA $RESCOL
    } | tee "$OUTPUT_FILE" >> "$INSTANCE_LOG"
    if [[ -f "$DIR_SALIDA/nuevos_camino.txt" ]]; then
      {
        echo "--- Nodos recorridos (semilla $CURRENT_SEMILLA) ---"
        cat "$DIR_SALIDA/nuevos_camino.txt"
        echo
      } >> "$INSTANCE_LOG"
    fi
    SEMILLA=$((SEMILLA + 1))
  done
done