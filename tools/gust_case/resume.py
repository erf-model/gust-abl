"""Which runs of a written case already hold their answer, so run.sh and the set-up window skip them.

A wind direction (one ``run`` of the builder's Result.runs: its folder, the column run and the patches) is
finished when

* its record ``.gust_ran`` holds the decks the last run of that direction used, the case's data list
  ``.gust_case_data`` (the SHA-256 of each data file the case tool wrote: terrain, roughness map, copied files)
  and the arguments it was given, and they equal the current ones (a case written again with other settings or
  other data files, or other ``key=value`` arguments to run.sh, runs again);
* the terrain run's log.txt and each patch's end in GUST's steady verdict (``gust: steady after N iterations``, or
  the joint run's ``gust: joint run: <name> steady after N joint iterations``), and a plotfile ``<io.plot_file><step>`` written after that verdict holds its Header.

A run that stopped at ``steady.max_iterations`` (``not steady``), diverged, was stopped or crashed, or whose
plotfile is gone, runs again. The record is written when a direction starts, after its old logs are removed, so a
direction stopped part-way never counts as finished.

run.sh holds the same test in shell (builder._run_script); tools/test_gust_resume.py checks that the two agree.
"""
import hashlib
import os
import re

RECORD = ".gust_ran"
DATA = ".gust_case_data"
# every verdict GUST prints (App.cpp), the joint run's too; only "steady" counts
VERDICT = re.compile(r"^gust: (?:joint run: \S+ )?(steady|not steady|diverged|unsteady run ended) after \d+ "
                     r"(?:joint )?iterations")


def plot_prefix(deck_path):
    """io.plot_file of a GUST deck (GUST's default ``plt`` when the deck leaves it out)."""
    prefix = "plt"
    with open(deck_path, encoding="utf-8", errors="replace") as handle:   # a comment in another encoding
        for line in handle:
            parts = line.split("#")[0].split()
            if len(parts) >= 3 and parts[0] == "io.plot_file" and parts[1] == "=":
                prefix = parts[2]
    return prefix


def steady_plotfile(folder):
    """The plotfile written after the last steady verdict in folder/log.txt, if it exists; else None."""
    log = os.path.join(folder, "log.txt")
    deck = os.path.join(folder, "inputs")
    if not (os.path.isfile(log) and os.path.isfile(deck)):
        return None
    prefix = plot_prefix(deck)
    steady, found = False, None
    with open(log, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            m = VERDICT.match(line)
            if m:
                steady, found = m.group(1) == "steady", None
                continue
            words = line.split()
            if (steady and len(words) >= 3 and words[0] == "gust:" and words[1] == "wrote"
                    and words[2].startswith(prefix) and words[2][len(prefix):].isdigit()
                    and words[2][len(prefix):].isascii()):
                found = words[2]
    if not steady or found is None or not os.path.isfile(os.path.join(folder, found, "Header")):
        return None
    return found


def decks(run):
    """The deck files of a direction, relative to its folder with ``/`` (run.sh uses them on every system)."""
    files = ["inputs"]
    if run["column"]:
        files.append("column/inputs")
    files += ["nest_%s/inputs" % p for p in run["patches"]]
    return files


def record_files(run):
    """The files whose bytes the record holds, relative to the direction's folder: the decks, then the case's
    data list (in the case folder, one level up for a sweep)."""
    return decks(run) + [DATA if run["dir"] == "." else "../" + DATA]


def write_data(case_dir, files):
    """Write the case's data list: the SHA-256 and name of each of ``files`` (relative to ``case_dir``)."""
    lines = ["# SHA-256 of the data files make_gust_case.py wrote (a changed file runs every direction again)\n"]
    for rel in files:
        digest = hashlib.sha256()
        with open(os.path.join(case_dir, rel), "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        lines.append("%s  %s\n" % (digest.hexdigest(), rel))
    with open(os.path.join(case_dir, DATA), "w", encoding="utf-8", newline="\n") as out:
        out.writelines(lines)


def record_bytes(base, run, args):
    """What the record of the direction in ``base`` holds for the current decks, data list and run.sh arguments
    ``args``: each file's bytes, then ``# arguments:`` and each argument after a space (run.sh writes the same
    bytes)."""
    out = b""
    for rel in record_files(run):
        with open(os.path.join(base, rel), "rb") as handle:
            out += handle.read()
    return out + b"# arguments:" + b"".join(b" " + os.fsencode(a) for a in args) + b"\n"


def finished(base, run, args=()):
    """True when the direction in ``base`` holds its answer for the current decks and ``args`` (module text)."""
    record = os.path.join(base, RECORD)
    if not os.path.isfile(record):
        return False
    try:
        with open(record, "rb") as handle:
            if handle.read() != record_bytes(base, run, args):
                return False
    except OSError:                                # a deck or the data list is missing: write the case again
        return False
    folders = [base] + [os.path.join(base, "nest_" + p) for p in run["patches"]]
    return all(steady_plotfile(f) is not None for f in folders)


def start(base, run, args=()):
    """A direction starts: remove its runs' old logs, then write its record."""
    for folder in [base] + [os.path.join(base, "nest_" + p) for p in run["patches"]]:
        log = os.path.join(folder, "log.txt")
        if os.path.isfile(log):
            os.remove(log)
    data = record_bytes(base, run, args)
    with open(os.path.join(base, RECORD), "wb") as handle:
        handle.write(data)
