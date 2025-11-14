import argparse
import math
import re
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import cm

Camino = Sequence[int]
MapaVisitas = Mapping[Tuple[int, int], int]


PROYECTO_RAIZ = Path(__file__).resolve().parent.parent
NOMBRE_ARCHIVO_POSICIONES = "posiciones.json"


def _resolver_entrada(ruta: Path) -> Path:
    """Resuelve rutas de entrada tolerando ejecuciones fuera de la raíz del repo."""

    if ruta.is_absolute():
        return ruta

    ruta_desde_cwd = Path.cwd() / ruta
    if ruta_desde_cwd.exists():
        return ruta_desde_cwd

    ruta_desde_raiz = PROYECTO_RAIZ / ruta
    if ruta_desde_raiz.exists():
        return ruta_desde_raiz

    return ruta_desde_cwd


def _resolver_salida(ruta: Path) -> Path:
    """Resuelve rutas de salida generando directorios relativos a la raíz si es necesario."""

    if ruta.is_absolute():
        return ruta

    ruta_desde_cwd = Path.cwd() / ruta
    if ruta_desde_cwd.exists():
        return ruta_desde_cwd

    return PROYECTO_RAIZ / ruta


def _ruta_posiciones(carpeta_salida: Path) -> Path:
    """Calcula la ruta donde se cachearán las posiciones de los nodos."""

    return carpeta_salida / NOMBRE_ARCHIVO_POSICIONES


def _guardar_posiciones(posiciones: Mapping[int, Tuple[float, float]], destino: Path) -> None:
    """Serializa las posiciones calculadas para reutilizarlas en ejecuciones futuras."""

    destino.parent.mkdir(parents=True, exist_ok=True)
    datos = {str(nodo): [float(x), float(y)] for nodo, (x, y) in posiciones.items()}
    destino.write_text(json.dumps(datos, indent=2), encoding="utf-8")


def _cargar_posiciones(ruta: Path) -> Optional[Dict[int, Tuple[float, float]]]:
    """Intenta cargar un diccionario de posiciones desde disco."""

    if not ruta.is_file():
        return None

    try:
        datos = json.loads(ruta.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None

    posiciones: Dict[int, Tuple[float, float]] = {}
    for clave, valor in datos.items():
        try:
            nodo = int(clave)
            x, y = float(valor[0]), float(valor[1])
        except (ValueError, TypeError, IndexError):
            return None
        posiciones[nodo] = (x, y)

    return posiciones


def _calcular_posiciones(
    grafo: nx.DiGraph, ruta_cache: Path, semilla: int = 42
) -> Dict[int, Tuple[float, float]]:
    """Obtiene las posiciones para los nodos, reutilizando un cache si está disponible."""

    posiciones = _cargar_posiciones(ruta_cache)
    if posiciones is not None:
        nodos_grafo = set(grafo.nodes)
        nodos_posiciones = set(posiciones.keys())
        if nodos_grafo == nodos_posiciones:
            return posiciones

    posiciones = nx.spring_layout(grafo, seed=semilla)
    _guardar_posiciones(posiciones, ruta_cache)
    return posiciones


def leer_aristas_req(ruta_instancia: Path) -> List[Tuple[str, int, int, float, float]]:
    """Lee las aristas requeridas desde una instancia y las retorna en bruto."""

    if not ruta_instancia.is_file():
        raise FileNotFoundError(f"No se encontró la instancia en {ruta_instancia}.")

    encabezado: Dict[str, str] = {}
    aristas: List[Tuple[str, int, int, float, float]] = []

    with ruta_instancia.open("r", encoding="utf-8") as archivo:
        for _ in range(8):
            linea = archivo.readline().strip()
            if not linea:
                continue
            clave, valor = linea.split(":")
            encabezado[clave.strip()] = valor.strip()

        archivo.readline()  # encabezado de sección "ARISTAS_REQ"

        cantidad_aristas = int(encabezado.get("ARISTAS_REQ", 0))
        for _ in range(cantidad_aristas):
            linea = archivo.readline().strip()
            if not linea:
                continue
            tipo, nodo1, nodo2, costo_recorrer, costo_recolectar = linea.split()
            aristas.append(
                (
                    tipo,
                    int(nodo1),
                    int(nodo2),
                    float(costo_recorrer),
                    float(costo_recolectar),
                )
            )

    return aristas


def construir_grafo_dirigido(aristas: Iterable[Tuple[str, int, int, float, float]]) -> nx.DiGraph:
    """Construye un grafo dirigido a partir de la lista de aristas de la instancia."""

    grafo = nx.DiGraph()

    for tipo, origen, destino, _costo_recorrido, _costo_recoleccion in aristas:
        grafo.add_edge(origen, destino)
        if tipo.lower() == "bi":
            grafo.add_edge(destino, origen)

    return grafo


def _extender_grafo_con_rutas(
    grafo: nx.DiGraph, rutas: Mapping[int, Camino]
) -> Tuple[Sequence[int], Sequence[Tuple[int, int]]]:
    """Agrega al grafo los nodos y aristas que aparecen en las rutas pero no en la instancia."""

    nodos_agregados: List[int] = []
    aristas_agregadas: List[Tuple[int, int]] = []

    for ruta in rutas.values():
        for nodo in ruta:
            if nodo not in grafo:
                grafo.add_node(nodo)
                nodos_agregados.append(nodo)

        for origen, destino in zip(ruta[:-1], ruta[1:]):
            if not grafo.has_edge(origen, destino):
                grafo.add_edge(origen, destino)
                aristas_agregadas.append((origen, destino))

    return nodos_agregados, aristas_agregadas


def leer_mapa_visitas(ruta_mapa: Path) -> MapaVisitas:
    """Lee ``nuevos_mapa.txt`` devolviendo un diccionario {(origen, destino): visitas}."""

    if not ruta_mapa.is_file():
        raise FileNotFoundError(f"No se encontró el mapa de visitas en {ruta_mapa}.")

    patron_entero = re.compile(r"-?\d+")
    visitas: Dict[Tuple[int, int], int] = {}

    with ruta_mapa.open("r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue

            if ":" not in linea:
                raise ValueError(
                    "Formato inválido en nuevos_mapa.txt; se esperaba 'origen,destino: cantidad'."
                )

            lado_arco, lado_valor = linea.split(":", maxsplit=1)
            nodos = patron_entero.findall(lado_arco)
            valores = patron_entero.findall(lado_valor)

            if len(nodos) < 2 or not valores:
                raise ValueError(
                    "Formato inválido en nuevos_mapa.txt; no se pudieron extraer nodos o cantidad."
                )

            origen, destino = int(nodos[0]), int(nodos[1])
            cantidad = int(valores[0])
            visitas[(origen, destino)] = cantidad

    return visitas


def leer_rutas_camiones(ruta_camino: Path) -> Dict[int, List[int]]:
    """Parsea ``nuevos_camino.txt`` agrupando los nodos por camión."""

    if not ruta_camino.is_file():
        raise FileNotFoundError(f"No se encontró el archivo de rutas en {ruta_camino}.")

    patron_entero = re.compile(r"-?\d+")
    rutas: Dict[int, List[int]] = {}

    camion_actual: Optional[int] = None
    buffer: List[str] = []

    def _guardar_camion() -> None:
        nonlocal camion_actual, buffer
        if camion_actual is None:
            return
        nodos = patron_entero.findall(" ".join(buffer))
        rutas[camion_actual] = [int(nodo) for nodo in nodos]
        buffer = []

    with ruta_camino.open("r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue

            if linea.lower().startswith("camion"):
                _guardar_camion()
                identificadores = patron_entero.findall(linea.split(":", maxsplit=1)[0])
                if not identificadores:
                    raise ValueError(
                        f"No se pudo determinar el número de camión en la línea: '{linea}'."
                    )
                camion_actual = int(identificadores[0])
                buffer.append(linea.split(":", maxsplit=1)[1])
            else:
                buffer.append(linea)

    _guardar_camion()

    return rutas


def _dibujar_base(
    grafo: nx.DiGraph,
    pos: Mapping[int, Tuple[float, float]],
    visitas: MapaVisitas,
    ax: plt.Axes,
) -> None:
    """Dibuja los nodos y arcos generales del grafo como contexto."""

    ax.axis("off")

    nx.draw_networkx_nodes(grafo, pos, node_color="#f7f7f7", edgecolors="#333333", linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(grafo, pos, font_size=8, ax=ax)

    edges_con_visita = [arco for arco in grafo.edges if visitas.get(arco, 0) > 0]
    edges_sin_visita = [arco for arco in grafo.edges if visitas.get(arco, 0) <= 0]

    if edges_sin_visita:
        nx.draw_networkx_edges(
            grafo,
            pos,
            edgelist=edges_sin_visita,
            edge_color="#c7c7c7",
            width=0.6,
            alpha=0.15,
            arrows=False,
            ax=ax,
        )

    if edges_con_visita:
        widths = [1.0 + math.log1p(visitas.get(arco, 0)) for arco in edges_con_visita]
        nx.draw_networkx_edges(
            grafo,
            pos,
            edgelist=edges_con_visita,
            edge_color="#7f7f7f",
            width=widths,
            alpha=0.45,
            arrows=False,
            ax=ax,
        )


def _dibujar_ruta_camion(
    grafo: nx.DiGraph,
    pos: Mapping[int, Tuple[float, float]],
    ruta: Camino,
    camion_id: int,
    carpeta_salida: Path,
    visitas: MapaVisitas,
    mostrar: bool,
) -> None:
    if len(ruta) < 2:
        raise ValueError(f"La ruta del camión {camion_id} debe contener al menos dos nodos.")

    figura, ax = plt.subplots(figsize=(12, 9))
    ax.set_title(f"Recorrido del camión {camion_id}")

    _dibujar_base(grafo, pos, visitas, ax)

    edges_ruta = list(zip(ruta[:-1], ruta[1:]))
    cmap = cm.get_cmap("cool", len(edges_ruta))

    for indice, arco in enumerate(edges_ruta):
        nx.draw_networkx_edges(
            grafo,
            pos,
            edgelist=[arco],
            edge_color=[cmap(indice)],
            width=3.0,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=10,
            ax=ax,
        )

    inicio, fin = ruta[0], ruta[-1]
    ax.scatter(*pos[inicio], s=180, c="#ff6f69", edgecolors="#222", zorder=5, label="Inicio")
    ax.scatter(*pos[fin], s=180, c="#88d8b0", edgecolors="#222", zorder=5, label="Fin")

    ax.legend(loc="upper right")

    carpeta_salida.mkdir(parents=True, exist_ok=True)
    archivo_salida = carpeta_salida / f"camion_{camion_id:02d}.png"
    figura.tight_layout()
    figura.savefig(archivo_salida, dpi=200)

    if mostrar:
        plt.show()

    plt.close(figura)


def visualizar_camiones(
    ruta_instancia: Path,
    ruta_camino: Path,
    ruta_mapa: Path,
    carpeta_salida: Path,
    mostrar: bool,
) -> None:
    aristas = leer_aristas_req(ruta_instancia)
    grafo = construir_grafo_dirigido(aristas)
    rutas_camiones = leer_rutas_camiones(ruta_camino)
    visitas = leer_mapa_visitas(ruta_mapa)

    nodos_extra, aristas_extra = _extender_grafo_con_rutas(grafo, rutas_camiones)
    if nodos_extra or aristas_extra:
        mensaje = [
            "Advertencia: las rutas incluyen elementos que no estaban en la instancia cargada.",
        ]
        if nodos_extra:
            mensaje.append(
                f"  - Nodos agregados dinámicamente: {sorted(set(nodos_extra))}"
            )
        if aristas_extra:
            mensaje.append(
                f"  - Aristas agregadas dinámicamente: {sorted(set(aristas_extra))}"
            )
        print("\n".join(mensaje))

    if grafo.number_of_nodes() == 0:
        raise ValueError("El grafo construido está vacío; revise la instancia proporcionada.")

    ruta_posiciones = _ruta_posiciones(carpeta_salida)
    pos = _calcular_posiciones(grafo, ruta_posiciones)

    for camion_id in sorted(rutas_camiones):
        ruta = rutas_camiones[camion_id]
        _dibujar_ruta_camion(grafo, pos, ruta, camion_id, carpeta_salida, visitas, mostrar)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualiza los recorridos de los camiones registrados por una hormiga.",
    )
    parser.add_argument(
        "instancia",
        type=Path,
        help="Ruta al archivo de instancia (por ejemplo, Instancias/QN/Seg_QN.txt)",
    )
    parser.add_argument(
        "--camino",
        type=Path,
        default=Path("resultados/nuevos_camino.txt"),
        help="Ruta a nuevos_camino.txt con los recorridos de los camiones.",
    )
    parser.add_argument(
        "--mapa",
        type=Path,
        default=Path("resultados/nuevos_mapa.txt"),
        help="Ruta a nuevos_mapa.txt con las frecuencias de los arcos visitados.",
    )
    parser.add_argument(
        "--salida",
        type=Path,
        default=Path("visualizaciones_camiones"),
        help="Carpeta donde se guardarán las imágenes generadas.",
    )
    parser.add_argument(
        "--mostrar",
        action="store_true",
        help="Muestra las gráficas en pantalla además de guardarlas en disco.",
    )

    args = parser.parse_args()

    ruta_instancia = _resolver_entrada(args.instancia)
    ruta_camino = _resolver_entrada(args.camino)
    ruta_mapa = _resolver_entrada(args.mapa)
    carpeta_salida = _resolver_salida(args.salida)

    visualizar_camiones(ruta_instancia, ruta_camino, ruta_mapa, carpeta_salida, args.mostrar)


if __name__ == "__main__":
    main()