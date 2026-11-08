"""
diagnostics/find_ew_bug_sites.py
AP: 2026-09-20

Locate every copy of the C IV equivalent-width continuum indexing defect in the
repository, so each can be inspected by hand.

The defect is `continuum[i]` inside a loop whose index runs over the EW-masked
subset (`wave_r[EW]`), while `continuum` itself is full length. The correct form,
which the plotting call in the same functions already uses, is `continuum[EW]`.

Run from the Richards root:
    python "HST Paper/diagnostics/find_ew_bug_sites.py"
"""

import glob
import json
import os
import re

BUG = re.compile(r"continuum\[i\]")
OK = re.compile(r"continuum\[EW\]")
SEP = chr(92)  # backslash, kept out of literals for shell-safety


def norm(p):
    return p.replace(SEP, "/")


def scan():
    rows = []
    for path in glob.glob("**/*.ipynb", recursive=True) + glob.glob("**/*.py", recursive=True):
        try:
            if path.endswith(".ipynb"):
                nb = json.load(open(path, encoding="utf-8"))
                for ci, c in enumerate(nb.get("cells", [])):
                    if c.get("cell_type") != "code":
                        continue
                    src = "".join(c.get("source", []))
                    if not BUG.search(src):
                        continue
                    hits = [(j + 1, l.strip()) for j, l in enumerate(src.split("\n"))
                            if BUG.search(l)]
                    rows.append((norm(path), "cell index %d (In[%s])"
                                 % (ci, c.get("execution_count")), hits, bool(OK.search(src))))
            else:
                src = open(path, encoding="utf-8", errors="replace").read()
                if not BUG.search(src):
                    continue
                hits = [(j + 1, l.strip()) for j, l in enumerate(src.split("\n"))
                        if BUG.search(l)]
                rows.append((norm(path), "file", hits, bool(OK.search(src))))
        except Exception:
            pass
    return rows


def bucket(p):
    if p.startswith("HST/"):
        return 0, "TREVOR'S UNTOUCHED ARCHIVE (HST/)"
    if p.startswith("Trevor Code/"):
        return 1, "TREVOR CODE/ (his code with AP edits)"
    return 2, "HST PAPER/ (the current pipeline)"


def main():
    rows = scan()
    rows.sort(key=lambda r: (bucket(r[0])[0], 1 if "checkpoint" in r[0] else 0, r[0]))
    last = None
    for p, loc, hits, has_ok in rows:
        b, label = bucket(p)
        if b != last:
            print("\n" + "=" * 78)
            print(label)
            print("=" * 78)
            last = b
        note = "   <- same cell also uses correct continuum[EW] for the shading" if has_ok else ""
        ckpt = "   (.ipynb_checkpoints autosave copy)" if "checkpoint" in p else ""
        print("\n%s%s" % (p, ckpt))
        print("    %s   %d occurrence(s)%s" % (loc, len(hits), note))
        for ln, txt in hits:
            print("        line %-3d of that cell/file: %s" % (ln, txt[:95]))
    print("\n%d file(s) contain the pattern." % len(rows))


if __name__ == "__main__":
    main()
