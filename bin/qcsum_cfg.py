#!/usr/bin/env python3
# Resolve panel settings and run qcsum_metrics.pl.
#   qcsum_cfg.py CONFIG PANEL get KEY
#   qcsum_cfg.py CONFIG PANEL perl|print-perl SCRIPT PREFIX QCFOLDER
import os, shlex, sys, yaml

PANEL_MAP = {
    "Comprehensive Heme Panel": "heme_comp",
    "MPN Screen Panel (JAK2, CALR, MPL)": "heme_mpn",
    "JAK2": "heme_jak2",
    "Peripheral Blood Lymphoma Panel": "heme_pblp",
    "Tumor": "pcp_tumor",
    "Normal": "pcp_normal",
}
PERL_KEYS = ["pipeline_version", "platform", "fail_min_align_pct", "covered",
             "fail_min_roi_pct", "fail_min_avgcov", "fail_min_reads", "capture"]

def get_config(path, panel):
    panel = PANEL_MAP.get(panel, panel)
    with open(path) as f:
        y = yaml.safe_load(f)
    if not y:
        sys.exit(f"Failed to load {path}. Ensure it is a valid YAML file.")
    if panel not in y:
        sys.exit(f"Panel '{panel}' not found in {path}.")
    if panel.startswith("heme_") and panel != "heme_comp":
        d = dict(y.get("heme_comp", {}))
        ti = (y.get(panel) or {}).get("target_intervals")
        if ti:
            d["target_intervals"] = ti
    else:
        d = y.get(panel) or {}
    return {k: str(v) for k, v in d.items()}

def need(cfg, k):
    if k not in cfg:
        sys.exit(f"qcsum config is missing key: {k}")
    return cfg[k]

def main():
    a = sys.argv[1:]
    if len(a) < 3:
        sys.exit("usage: qcsum_cfg.py CONFIG PANEL get KEY | perl SCRIPT PREFIX QCFOLDER")
    cfg = get_config(a[0], a[1])
    if a[2] == "get" and len(a) == 4:
        print(need(cfg, a[3]))
    elif a[2] in {"perl", "print-perl"} and len(a) == 6:
        args = ["perl", a[3], "--prefix", a[4], "--qcfolder", a[5]]
        for k in PERL_KEYS:
            args += ["--" + k, need(cfg, k)]
        if a[2] == "print-perl":
            print(shlex.join(args))
        else:
            sys.stdout.flush()
            os.execvp("perl", args)
    else:
        sys.exit("bad arguments: " + " ".join(a))

main()
