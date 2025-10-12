#!/usr/bin/env python3
"""
Trova loop nelle dipendenze locali (#include "...") di un progetto C/C++.
Autore: ChatGPT (2025)
"""

import os
import re
import sys
from collections import defaultdict

try:
    import networkx as nx
except ImportError:
    print("❌ Devi installare networkx prima di eseguire questo script:")
    print("   pip install networkx")
    sys.exit(1)

# Pattern: match solo #include "..."
INCLUDE_PATTERN = re.compile(r'^\s*#include\s*"([^"]+)"')


def find_includes(root_dir):
    """Ritorna un dizionario file -> set(includes)"""
    deps = defaultdict(set)
    for dirpath, _, files in os.walk(root_dir):
        for f in files:
            if f.endswith((".h", ".hpp", ".tpp", ".cpp")):
                full = os.path.join(dirpath, f)
                rel = os.path.relpath(full, root_dir)
                try:
                    with open(full, "r", encoding="utf-8", errors="ignore") as fh:
                        for line in fh:
                            m = INCLUDE_PATTERN.match(line)
                            if m:
                                deps[rel].add(m.group(1))
                except Exception as e:
                    print(f"⚠️ Errore leggendo {full}: {e}")
    return deps


def build_graph(deps):
    """Costruisce un grafo diretto NetworkX a partire dal dizionario deps"""
    G = nx.DiGraph()
    for src, includes in deps.items():
        for inc in includes:
            G.add_edge(src, inc)
    return G


def find_loops(G):
    """Trova tutti i cicli semplici nel grafo"""
    return list(nx.simple_cycles(G))


def main():
    if len(sys.argv) < 2:
        print("Uso: python find_dependency_loops.py <path_progetto>")
        sys.exit(1)

    root = sys.argv[1]
    print(f"🔍 Analizzo dipendenze in: {root}")

    deps = find_includes(root)
    G = build_graph(deps)
    loops = find_loops(G)

    print(f"\n📦 File analizzati: {len(deps)}")
    print(f"🔁 Cicli trovati: {len(loops)}\n")

    if not loops:
        print("✅ Nessun loop trovato! Tutto pulito 🚀")
    else:
        for i, loop in enumerate(loops, start=1):
            chain = " → ".join(loop + [loop[0]])
            print(f"🔁 Loop #{i}: {chain}")


if __name__ == "__main__":
    main()
