
import argparse
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# ---------------------------------------------------------------------------
# Datos auxiliares

EdgeKey = Tuple[int, int]


@dataclass
class RunArtifacts:
    nombre: str
    directorio: Path
    archivo_camino: Path
    archivo_mapa: Optional[Path]


# ---------------------------------------------------------------------------
# Lectura de argumentos


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera imágenes de los recorridos construidos por RESCOL a partir "
            "de los archivos *_camino.txt y *_mapa.txt almacenados en un "
            "directorio de resultados."
        )
    )
    parser.add_argument(
        "--dir-resultados",
        type=Path,
        default=Path("./resultados"),
        help="Directorio donde se almacenan los resultados de RESCOL.",
    )
    parser.add_argument(
        "--prefijo",
        default="resultado",
        help=(
            "Prefijo de los subdirectorios que contienen corridas. Solo se "
            "procesarán carpetas cuyo nombre comience con este valor."
        ),
    )
    parser.add_argument(
        "--dir-salida",
        type=Path,
        default=Path("./imagenes_recorridos"),
        help="Carpeta donde se guardarán las imágenes generadas.",
    )
    parser.add_argument(
        "--layout",
        choices=("spring", "kamada_kawai", "planar"),
        default="kamada_kawai",
        help="Layout de NetworkX utilizado para posicionar los nodos.",
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=12.0,
        help="Ancho de la figura (en pulgadas).",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=8.0,
        help="Alto de la figura (en pulgadas).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolución de la imagen de salida (puntos por pulgada).",
    )
    parser.add_argument(
        "--mostrar",
        action="store_true",
        help="Muestra la figura en pantalla además de guardarla (requiere entorno gráfico).",
    )
    parser.add_argument(
        "--log-nivel",
        choices=("debug", "info", "warning", "error"),
        default="info",
        help="Nivel de detalle de los mensajes de log.",
    )

    args = parser.parse_args(argv)
    if args.fig_width <= 0 or args.fig_height <= 0:
        raise SystemExit("El ancho y alto de la figura deben ser positivos.")
    if args.dpi <= 0:
        raise SystemExit("El valor de dpi debe ser positivo.")
    return args


# ---------------------------------------------------------------------------
# Funciones utilitarias para localizar artefactos


def discover_runs(base_dir: Path, prefix: str) -> List[RunArtifacts]:
    if not base_dir.exists():
        raise FileNotFoundError(f"El directorio de resultados '{base_dir}' no existe.")

    runs: List[RunArtifacts] = []
    for child in sorted(base_dir.iterdir()):
        if not child.is_dir():
            continue
        if prefix and not child.name.startswith(prefix):
            logging.debug("Se ignora '%s' porque no coincide con el prefijo '%s'.", child, prefix)
            continue

        camino = locate_file(child, ["camino"])
        mapa = locate_file(child, ["mapa"])
        if camino is None:
            logging.warning(
                "No se encontró un archivo de camino en '%s'; se omite la corrida.",
                child,
            )
            continue
        runs.append(
            RunArtifacts(
                nombre=child.name,
                directorio=child,
                archivo_camino=camino,
                archivo_mapa=mapa,
            )
        )
    return runs


def locate_file(root: Path, keywords: Iterable[str]) -> Optional[Path]:
    keywords_normalized = [kw.lower() for kw in keywords]
    candidates: List[Path] = []
    for path in sorted(root.rglob("*.txt")):
        lower_name = path.name.lower()
        if all(kw in lower_name for kw in keywords_normalized):
            candidates.append(path)
    if not candidates:
        return None
    # Prefer archivos ubicados directamente en el directorio raíz
    candidates.sort(key=lambda p: (p.parent != root, p))
    return candidates[0]


# ---------------------------------------------------------------------------
# Lectura de archivos de resultados


def read_route(path: Path) -> List[int]:
    route: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            numbers = re.findall(r"-?\d+", line)
            if not numbers:
                logging.debug("No se encontraron números en la línea '%s' del archivo '%s'.", line, path)
                continue
            for token in numbers:
                route.append(int(token))
    if len(route) < 2:
        raise ValueError(
            f"El archivo '{path}' no contiene al menos dos nodos para definir un recorrido."
        )
    return route


def read_edge_counts(path: Path) -> Dict[EdgeKey, int]:
    if path is None:
        return {}
    edge_counts: Dict[EdgeKey, int] = {}
    pattern = re.compile(r"\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*:?\s*(-?\d+)")
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            match = pattern.search(line)
            if not match:
                logging.debug("Formato desconocido en '%s': '%s'", path, line)
                continue
            origen, destino, cantidad = map(int, match.groups())
            edge_counts[(origen, destino)] = cantidad
    return edge_counts


# ---------------------------------------------------------------------------
# Visualización


def create_layout(graph: nx.DiGraph, layout: str) -> Dict[int, np.ndarray]:
    if layout == "spring":
        return nx.spring_layout(graph, seed=42)
    if layout == "kamada_kawai":
        return nx.kamada_kawai_layout(graph)
    if layout == "planar":
        try:
            return nx.planar_layout(graph)
        except (nx.NetworkXException, ValueError):
            logging.warning("No fue posible aplicar el layout planar; se usará spring_layout.")
            return nx.spring_layout(graph, seed=42)
    raise ValueError(f"Layout desconocido: {layout}")


def build_graph(edge_counts: Dict[EdgeKey, int], route: Sequence[int]) -> nx.DiGraph:
    graph = nx.DiGraph()
    if edge_counts:
        for (origen, destino), cantidad in edge_counts.items():
            graph.add_edge(origen, destino, cantidad=cantidad)
    # Asegurar que todos los nodos del recorrido estén presentes
    for node in route:
        graph.add_node(node)
    if graph.number_of_edges() == 0:
        # Si no hay información de arcos, construimos un camino simple
        graph.add_edges_from(zip(route[:-1], route[1:]), cantidad=1)
    return graph


def draw_run(
    artifacts: RunArtifacts,
    layout_name: str,
    fig_size: Tuple[float, float],
    dpi: int,
    mostrar: bool,
    output_dir: Path,
) -> Path:
    route = read_route(artifacts.archivo_camino)
    edge_counts = read_edge_counts(artifacts.archivo_mapa) if artifacts.archivo_mapa else {}
    graph = build_graph(edge_counts, route)

    pos = create_layout(graph, layout_name)
    figure = plt.figure(figsize=fig_size)
    ax = figure.add_subplot(1, 1, 1)
    ax.set_axis_off()

    base_edges = list(graph.edges())
    path_edges = list(zip(route[:-1], route[1:]))

    # Dibuja primero todos los arcos en color tenue
    base_widths = [0.8 + math.log(graph[u][v].get("cantidad", 1) + 1) for u, v in base_edges]
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_color="#fee08b",
        edgecolors="#636363",
        linewidths=0.8,
        ax=ax,
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        edgelist=base_edges,
        width=base_widths,
        edge_color="#bdbdbd",
        arrows=False,
        alpha=0.5,
        ax=ax,
    )

    # Destaca el recorrido en un degradado de colores
    if path_edges:
        colores = cm.get_cmap("viridis", len(path_edges))
        for idx, edge in enumerate(path_edges):
            nx.draw_networkx_edges(
                graph,
                pos,
                edgelist=[edge],
                width=2.5,
                edge_color=[colores(idx)],
                arrows=True,
                arrowstyle="-|>",
                arrowsize=12,
                ax=ax,
            )

    # Etiquetas de nodos
    nx.draw_networkx_labels(graph, pos, font_size=9, font_color="#1a1a1a", ax=ax)

    # Etiquetas con la cantidad de visitas por arco si están disponibles
    if edge_counts:
        etiquetas = {(u, v): str(edge_counts.get((u, v), "")) for u, v in base_edges}
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=etiquetas, font_size=7, ax=ax)

    ax.set_title(f"Recorrido: {artifacts.nombre}")
    plt.tight_layout()

    imagen_path = output_dir / f"{artifacts.nombre}_recorrido.png"
    plt.savefig(imagen_path, dpi=dpi, bbox_inches="tight")
    if mostrar:
        plt.show()
    plt.close(figure)
    return imagen_path


# ---------------------------------------------------------------------------
# Flujo principal


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_nivel.upper()),
        format="[%(levelname)s] %(message)s",
    )

    args.dir_salida.mkdir(parents=True, exist_ok=True)
    runs = discover_runs(args.dir_resultados, args.prefijo)
    if not runs:
        logging.warning(
            "No se encontraron corridas en '%s' con el prefijo '%s'.", args.dir_resultados, args.prefijo
        )
        return

    fig_size = (args.fig_width, args.fig_height)
    for artifacts in runs:
        logging.info("Procesando corrida '%s'", artifacts.nombre)
        try:
            image_path = draw_run(
                artifacts,
                layout_name=args.layout,
                fig_size=fig_size,
                dpi=args.dpi,
                mostrar=args.mostrar,
                output_dir=args.dir_salida,
            )
        except Exception as exc:  # pylint: disable=broad-except
            logging.error("No fue posible generar la imagen para '%s': %s", artifacts.nombre, exc)
            continue
        logging.info("Imagen guardada en '%s'", image_path)


if __name__ == "__main__":
    main()