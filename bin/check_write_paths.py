#!/usr/bin/env python3
# Defence in depth for --forbid_prefixes: refuse (exit 1) to write a path that is, or resolves through a symlink to,
# a location equal to or inside a forbidden prefix. Called at the top of every process that writes, before any write.
#   check_write_paths.py FORBID_CSV [--tree DIR] PATH...
# PATH may not exist yet; its existing part is resolved (os.path.realpath). --tree DIR also checks every entry below
# DIR (symlinks are checked, not followed). An empty FORBID_CSV checks nothing (production default).
import os, sys

def forms(p):
    a = os.path.abspath(p)
    return {a, os.path.realpath(a)}

def hit(path, prefixes):
    for x in prefixes:
        for y in forms(x):
            y = y.rstrip("/") or "/"
            for q in forms(path):
                if y == "/" or q == y or q.startswith(y + "/"):
                    return x, q
    return None

def main(a):
    if not a:
        sys.exit("usage: check_write_paths.py FORBID_CSV [--tree DIR] PATH...")
    forbid = [x.strip() for x in a[0].split(",") if x.strip()]
    paths, i = [], 1
    while i < len(a):
        if a[i] == "--tree" and i + 1 < len(a):
            root = a[i + 1]
            paths.append(root)
            if os.path.isdir(root) and not os.path.islink(root):
                for dp, ds, fs in os.walk(root):
                    paths += [os.path.join(dp, n) for n in ds + fs]
            i += 2
        else:
            paths.append(a[i])
            i += 1
    if not forbid:
        return 0
    bad = 0
    for p in paths:
        h = hit(p, forbid)
        if h:
            bad += 1
            print(f"ERROR: refusing to write {p}: it resolves to {h[1]}, equal to or inside forbidden prefix {h[0]} "
                  f"(--forbid_prefixes {a[0]}){' [symlink]' if os.path.islink(p) else ''}", file=sys.stderr)
    return 1 if bad else 0

sys.exit(main(sys.argv[1:]))
