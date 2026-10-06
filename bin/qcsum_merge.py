#!/usr/bin/env python3
"""Merge sample QC tables in sampleinfo order, keeping a single header."""
import sys
from pathlib import Path

output, *inputs = sys.argv[1:]
if not inputs:
    sys.exit('No sample QC tables to merge')
tables = [Path(path).read_text().splitlines() for path in inputs]
header = tables[0][0] if tables[0] else None
if any(len(table) != 2 or table[0] != header for table in tables):
    sys.exit('Invalid or inconsistent sample QC tables')
Path(output).write_text('\n'.join([header, *(table[1] for table in tables), '']) + '\n')
