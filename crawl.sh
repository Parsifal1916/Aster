#!/usr/bin/env python3
"""
Stampa un albero delle dipendenze (#include "...") nel terminale.
"""
import os
import re
from collections import defaultdict

pattern = re.compile(r'^\s*#include\s*"([^"]+)"')

def find_includes(root_dir):
    deps = defaultdict(set)
    for dirpath, _, files in os.walk(root_dir):
        for f in files:
            if f.endswith((".h", ".hpp", ".tpp", ".cpp")):
                full = os.path.join(dirpath, f)
                rel = os.path.relpath(full, root_dir)
                with open(full, 'r', encoding='utf-8', errors='ignore') as fh:
                    for line in fh:
                        m = pattern.match(line)
                        if m:
                            deps[rel].add(m.group(1))
    return deps


def print_tree(node, deps, visited=None, prefix=""):
    if visited is None:
        visited = set()
    if node in visited:
        print(prefix + f"↩ {node} (ciclo)")
        return
    visited.add(node)
    print(prefix + f"📄 {node}")
    children = sorted(deps.get(node, []))
    for i, child in enumerate(children):
        is_last = (i == len(children) - 1)
        connector = "└── " if is_last else "├── "
        subprefix = "    " if is_last else "│   "
        print_tree(child, deps, visited.copy(), prefix + connector)


def main():
    import sys
    root_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    deps = find_includes(root_dir)
    print(f"📦 Dipendenze in {root_dir}\n")
    for node in sorted(deps.keys()):
        print_tree(node, deps)
        print("")


if __name__ == "__main__":
    main()
