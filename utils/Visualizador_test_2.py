#from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

# Tipos auxiliares
Nodo = int
Arco = Tuple[Nodo, Nodo]


def parsear_camino(path: Path) -> List[Tuple[int, List[Nodo]]]:
    """Devuelve una lista ordenada de rutas por camión.

    Cada elemento de la lista corresponde a un camión e incluye su identificador
    (si está presente en el archivo) y la secuencia de nodos visitados.
    """

    rutas: List[Tuple[int, List[Nodo]]] = []
    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo con los caminos: {path}")

    patron_camion = re.compile(r"camion\s*(\d+)", re.IGNORECASE)

    with path.open("r", encoding="utf-8") as archivo:
        for indice, linea in enumerate(archivo, start=1):
            contenido = linea.strip()
            if not contenido:
                continue

            identificador = indice
            coincidencia = patron_camion.search(contenido)
            if coincidencia:
                identificador = int(coincidencia.group(1))
                _, _, contenido = contenido.partition(":")

            numeros = re.findall(r"\d+", contenido)
            ruta = [int(valor) for valor in numeros]
            if not ruta:
                # Si no hay nodos, ignoramos la línea para evitar generar imágenes vacías.
                continue
            rutas.append((identificador, ruta))

    if not rutas:
        raise ValueError(
            f"No se encontraron recorridos válidos en el archivo {path}. "
            "Compruebe que el algoritmo haya generado rutas para los camiones."
        )

    return rutas


def parsear_mapa(path: Path) -> Dict[Arco, int]:
    """Lee las repeticiones de cada arco visitado."""

    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo del mapa: {path}")

    conteos: Dict[Arco, int] = {}
    patron_linea = re.compile(r"\((\d+)\s*,\s*(\d+)\)\s*:\s*(\d+)")

    with path.open("r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue
            coincidencia = patron_linea.search(linea)
            if not coincidencia:
                continue
            origen, destino, cantidad = map(int, coincidencia.groups())
            conteos[(origen, destino)] = cantidad

    if not conteos:
        raise ValueError(
            f"El archivo {path} no contiene arcos con recorridos. "
            "Ejecute el algoritmo con la opción --salida-dijkstra para generar datos."
        )

    return conteos


def construir_grafo(rutas: Sequence[Tuple[int, Sequence[Nodo]]], conteos: Dict[Arco, int]) -> nx.DiGraph:
    """Construye un grafo dirigido con la información de rutas y conteos de arcos."""

    grafo = nx.DiGraph()

    # Añadimos arcos conocidos a partir del mapa.
    for (origen, destino), cantidad in conteos.items():
        grafo.add_edge(origen, destino, count=cantidad)

    # Aseguramos que todos los nodos y arcos de las rutas estén presentes.
    for _, ruta in rutas:
        for nodo in ruta:
            grafo.add_node(nodo)
        for origen, destino in zip(ruta, ruta[1:]):
            atributos = grafo.get_edge_data(origen, destino, default={})
            if "count" not in atributos:
                atributos["count"] = conteos.get((origen, destino), 0)
            grafo.add_edge(origen, destino, **atributos)

    return grafo


def _obtener_posiciones(grafo: nx.DiGraph, layout: str) -> Dict[Nodo, Tuple[float, float]]:
    """Calcula posiciones estables para el grafo."""

    if layout == "kamada_kawai":
        return nx.kamada_kawai_layout(grafo)
    if layout == "planar":
        try:
            return nx.planar_layout(grafo)
        except nx.NetworkXException:
            # Si no es planar, se usa spring layout.
            pass
    # Valor por defecto y método de respaldo.
    return nx.spring_layout(grafo, seed=42)


def _normalizar_anchos(valores: Iterable[int], minimo: float = 1.0, maximo: float = 6.0) -> Dict[int, float]:
    """Normaliza los valores de conteo a un rango de anchos de arista."""

    valores = list(valores)
    if not valores:
        return {}
    minimo_valor = min(valores)
    maximo_valor = max(valores)
    if minimo_valor == maximo_valor:
        return {valor: (minimo + maximo) / 2.0 for valor in valores}

    escala = maximo - minimo
    return {valor: minimo + ((valor - minimo_valor) / escala) * (maximo - minimo) for valor in valores}


def dibujar_mapa_general(
    grafo: nx.DiGraph,
    posiciones: Dict[Nodo, Tuple[float, float]],
    conteos: Dict[Arco, int],
    destino_archivo: Path,
) -> None:
    """Genera una imagen PNG con el mapa completo de recorridos."""

    destino_archivo.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    ax.set_title("Mapa de arcos recorridos")
    ax.axis("off")

    nx.draw_networkx_nodes(grafo, posiciones, node_size=700, node_color="#f8f9fa", edgecolors="#343a40", linewidths=1.5)

    conteos_existentes = [grafo[u][v]["count"] for u, v in grafo.edges()]
    anchos = _normalizar_anchos(conteos_existentes)
    colores = [
        "#495057"
        if conteo > 0
        else "#ced4da"
        for conteo in conteos_existentes
    ]
    nx.draw_networkx_edges(
        grafo,
        posiciones,
        edgelist=list(grafo.edges()),
        width=[anchos[conteo] for conteo in conteos_existentes],
        edge_color=colores,
        arrows=True,
        arrowsize=15,
        arrowstyle="-|>",
        connectionstyle="arc3,rad=0.08",
    )

    etiquetas = {arco: str(int(grafo[arco[0]][arco[1]]["count"])) for arco in grafo.edges() if grafo[arco[0]][arco[1]]["count"] > 0}
    nx.draw_networkx_labels(grafo, posiciones, font_size=10, font_color="#212529")
    nx.draw_networkx_edge_labels(grafo, posiciones, edge_labels=etiquetas, font_color="#1c7ed6")

    plt.tight_layout()
    plt.savefig(destino_archivo, dpi=200)
    plt.close()


def dibujar_camion(
    grafo: nx.DiGraph,
    posiciones: Dict[Nodo, Tuple[float, float]],
    conteos: Dict[Arco, int],
    identificador: int,
    ruta: Sequence[Nodo],
    destino_archivo: Path,
) -> None:
    """Genera la imagen de un camión resaltando su recorrido."""

    destino_archivo.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    ax.set_title(f"Camión {identificador}")
    ax.axis("off")

    # Base del grafo
    nx.draw_networkx_nodes(grafo, posiciones, node_size=650, node_color="#ffffff", edgecolors="#343a40", linewidths=1.5)
    nx.draw_networkx_edges(
        grafo,
        posiciones,
        edgelist=list(grafo.edges()),
        width=1,
        alpha=0.2,
        edge_color="#868e96",
        arrows=False,
    )
    nx.draw_networkx_labels(grafo, posiciones, font_size=10, font_color="#212529")

    arcos_ruta = list(zip(ruta, ruta[1:]))
    if not arcos_ruta:
        plt.savefig(destino_archivo, dpi=200)
        plt.close()
        return

    max_conteo = max(conteos.values() or [1])
    anchos = []
    for arco in arcos_ruta:
        conteo = conteos.get(arco, 0)
        if max_conteo == 0:
            ancho = 3.0
        else:
            ancho = 2.0 + 4.0 * (conteo / max_conteo)
        anchos.append(ancho)

    nx.draw_networkx_edges(
        grafo,
        posiciones,
        edgelist=arcos_ruta,
        width=anchos,
        edge_color="#d9480f",
        arrows=True,
        arrowsize=22,
        arrowstyle="-|>",
        connectionstyle="arc3,rad=0.1",
    )

    # Etiqueta de orden del recorrido
    for indice, (origen, destino) in enumerate(arcos_ruta, start=1):
        x = (posiciones[origen][0] + posiciones[destino][0]) / 2
        y = (posiciones[origen][1] + posiciones[destino][1]) / 2
        ax.text(x, y, str(indice), color="#212529", fontsize=9, bbox=dict(boxstyle="round,pad=0.2", facecolor="#fff3bf", alpha=0.8))

    plt.tight_layout()
    plt.savefig(destino_archivo, dpi=200)
    plt.close()


def generar_imagenes(
    directorio_resultados: Path,
    prefijo: str,
    directorio_salida: Path,
    layout: str = "spring",
) -> None:
    """Genera todas las imágenes necesarias a partir de los resultados."""

    archivo_camino = directorio_resultados / f"{prefijo}_camino.txt"
    archivo_mapa = directorio_resultados / f"{prefijo}_mapa.txt"

    rutas = parsear_camino(archivo_camino)
    conteos = parsear_mapa(archivo_mapa)

    grafo = construir_grafo(rutas, conteos)
    posiciones = _obtener_posiciones(grafo, layout)

    destino_mapa = directorio_salida / f"{prefijo}_mapa_general.png"
    dibujar_mapa_general(grafo, posiciones, conteos, destino_mapa)

    for identificador, ruta in rutas:
        destino = directorio_salida / f"{prefijo}_camion_{identificador}.png"
        dibujar_camion(grafo, posiciones, conteos, identificador, ruta, destino)


def crear_argumentos() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Genera imágenes PNG para los recorridos de camiones a partir de los "
            "archivos <prefijo>_camino.txt y <prefijo>_mapa.txt producidos por RESCOL."
        )
    )
    parser.add_argument(
        "--dir-resultados",
        type=Path,
        default=Path("./resultados"),
        help="Directorio que contiene los archivos de salida del algoritmo.",
    )
    parser.add_argument(
        "--prefijo",
        type=str,
        default="nuevos",
        help="Prefijo utilizado al generar los archivos de salida.",
    )
    parser.add_argument(
        "--dir-salida",
        type=Path,
        default=Path("./imagenes_recorridos"),
        help="Directorio donde se guardarán las imágenes generadas.",
    )
    parser.add_argument(
        "--layout",
        choices=["spring", "kamada_kawai", "planar"],
        default="spring",
        help="Método de distribución para posicionar los nodos en el plano.",
    )
    return parser


def main() -> None:
    parser = crear_argumentos()
    argumentos = parser.parse_args()
    generar_imagenes(argumentos.dir_resultados, argumentos.prefijo, argumentos.dir_salida, argumentos.layout)


if __name__ == "__main__":
    main()