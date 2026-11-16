"""Tests of the GUST input checker (tools/check_gust_inputs.py, tools/gust_input_checker).

    cd tools && python3 -m unittest -v test_gust_input_checker

Three parts:

1. the checker's own self-test table (gust_input_checker/selftest.py): every
   check has a deck that trips it and a sibling that does not;
2. the decks the CTests run: every test of Tests/CMakeLists.txt that runs gust
   and must start is checked with its command-line keys, and must have no
   ERROR (missing files that the test creates at run time excepted);
3. the decks the CTests expect to stop (the abort_* tests): the checker must
   report an ERROR on a key the expected message names.

Tests/CMakeLists.txt is read by a small CMake reader (set, foreach, if) that
resolves the variables the tests use; a test it cannot resolve is counted as
skipped, and the counts have a floor so that the sweeps cannot pass empty.
"""
import os
import re
import unittest

import gust_input_checker as gic
from gust_input_checker import Options, Severity, selftest
from gust_input_checker.checks import values
from gust_input_checker.context import OWNED
from gust_input_checker.reference import load_reference

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TESTS = os.path.join(REPO, "Tests")
GUST = "$<TARGET_FILE:gust>"
# PREFIX scripts that run the command as CMake wrote it (read: each calls subprocess.run(cmd) with no added key).
TRANSPARENT = ("check_mesh.py", "check_steady_hill.py", "check_steady_column.py", "check_column_target.py",
               "check_steady_neutral.py")
BUILD = "<build>"


# ---------------------------------------------------------------------------
# A reader of the part of CMake that Tests/CMakeLists.txt uses
# ---------------------------------------------------------------------------
def cmake_commands(text):
    """[(name, [(text, quoted)]), ...]: the commands of a CMake file, comments removed."""
    out, i, n = [], 0, len(text)
    ident = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
    while i < n:
        c = text[i]
        if c == "#":
            i = text.find("\n", i)
            i = n if i < 0 else i
            continue
        m = ident.match(text, i)
        if not m:
            i += 1
            continue
        name, j = m.group(0), m.end()
        while j < n and text[j] in " \t":
            j += 1
        if j >= n or text[j] != "(":
            i = m.end()
            continue
        args, depth, j, cur = [], 1, j + 1, ""
        while j < n and depth:
            c = text[j]
            if c == '"':
                k, s = j + 1, ""
                while k < n and text[k] != '"':
                    if text[k] == "\\" and k + 1 < n:
                        s += {"n": "\n", "t": "\t"}.get(text[k + 1], text[k + 1])
                        k += 2
                    else:
                        s += text[k]
                        k += 1
                if cur:
                    args.append((cur, False))
                    cur = ""
                args.append((s, True))
                j = k + 1
                continue
            if c == "#":
                j = text.find("\n", j)
                j = n if j < 0 else j
                continue
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            if c in " \t\r\n" and depth == 1:
                if cur:
                    args.append((cur, False))
                    cur = ""
            else:
                cur += c
            j += 1
        if cur:
            args.append((cur, False))
        out.append((name.lower(), args))
        i = j + 1
    return out


class CMakeTests:
    """Runs the commands and keeps the tests: {name: {"prefix": [...], "command": [...]}} in file order."""

    def __init__(self, path):
        self.vars = {"CMAKE_CURRENT_SOURCE_DIR": [os.path.dirname(path)], "PROJECT_SOURCE_DIR": [REPO],
                     "CMAKE_CURRENT_BINARY_DIR": [BUILD + "/Tests"], "CMAKE_BINARY_DIR": [BUILD],
                     "Python3_EXECUTABLE": ["python3"], "GUST_MPI": ["ON"], "GUST_PRECISION": ["DOUBLE"],
                     "GUST_GPU_BACKEND": ["NONE"], "CMAKE_COMMAND": ["cmake"]}
        self.tests = {}
        with open(path) as fh:
            self.run(cmake_commands(fh.read()))

    def expand(self, text):
        for _ in range(20):
            m = re.search(r"\$\{([A-Za-z0-9_]+)\}", text)
            if not m or m.group(1) not in self.vars:
                break
            text = text[:m.start()] + ";".join(self.vars[m.group(1)]) + text[m.end():]
        return text

    def args(self, raw):
        out = []
        for text, quoted in raw:
            text = self.expand(text)
            if quoted:
                out.append(text)
            else:
                out.extend(t for t in text.split(";") if t)
        return out

    def truth(self, raw):
        a = self.args(raw)

        def value(tok):
            return ";".join(self.vars[tok]) if tok in self.vars else tok

        def term(k):
            if a[k] == "NOT":
                v, k = term(k + 1)
                return (not v), k
            if a[k] == "TARGET":
                return False, k + 2
            if a[k] == "EXISTS":
                return os.path.exists(a[k + 1]), k + 2
            if k + 2 < len(a) + 0 and a[k + 1] == "STREQUAL":
                return value(a[k]) == value(a[k + 2]), k + 3
            v = value(a[k]) if a[k] in self.vars else ""
            return v.upper() not in ("", "0", "OFF", "FALSE", "NO", "NOTFOUND", "NONE") and True, k + 1

        v, k = term(0)
        while k < len(a):
            w, k2 = term(k + 1)
            v = (v and w) if a[k] == "AND" else (v or w)
            k = k2
        return v

    @staticmethod
    def block_end(cmds, i, opener, closer):
        depth = 0
        for j in range(i, len(cmds)):
            if cmds[j][0] == opener:
                depth += 1
            elif cmds[j][0] == closer:
                depth -= 1
                if depth == 0:
                    return j
        raise ValueError("no %s for the %s at command %d" % (closer, opener, i))

    def run(self, cmds):
        i = 0
        while i < len(cmds):
            name, raw = cmds[i]
            if name in ("function", "macro"):
                i = self.block_end(cmds, i, name, "end" + name) + 1
            elif name == "foreach":
                end = self.block_end(cmds, i, "foreach", "endforeach")
                a = self.args(raw)
                items = a[1:]
                if len(a) > 2 and a[1] == "IN" and a[2] == "LISTS":
                    items = [v for lst in a[3:] for v in self.vars.get(lst, [])]
                for item in items:
                    self.vars[a[0]] = [item]
                    self.run(cmds[i + 1:end])
                i = end + 1
            elif name == "if":
                end = self.block_end(cmds, i, "if", "endif")
                branches, depth, start, cond = [], 0, i + 1, raw
                for j in range(i + 1, end + 1):
                    if cmds[j][0] == "if":
                        depth += 1
                    elif cmds[j][0] == "endif" and depth:
                        depth -= 1
                    elif depth == 0 and cmds[j][0] in ("elseif", "else", "endif"):
                        branches.append((cond, start, j))
                        cond, start = (cmds[j][1] if cmds[j][0] == "elseif" else None), j + 1
                for cond, lo, hi in branches:
                    if cond is None or self.truth(cond):
                        self.run(cmds[lo:hi])
                        break
                i = end + 1
            else:
                if name == "set" and raw:
                    a = self.args(raw)
                    vals = a[1:]
                    if "CACHE" in vals:
                        vals = vals[:vals.index("CACHE")]
                    self.vars[a[0]] = vals
                elif name == "list" and len(raw) > 2 and raw[0][0] == "APPEND":
                    a = self.args(raw)
                    self.vars.setdefault(a[1], []).extend(a[2:])
                elif name in ("gust_add_test", "add_test"):
                    self.add(name, self.args(raw))
                i += 1

    def add(self, kind, a):
        words = ("NAME", "NP", "LABELS", "COMMAND", "PREFIX", "FIXTURES_SETUP", "FIXTURES_REQUIRED", "RUN_DIR",
                 "WORKING_DIRECTORY")
        if kind == "add_test":
            words = ("NAME", "COMMAND", "WORKING_DIRECTORY")
        parts, cur = {}, None
        for tok in a:
            if tok in words:
                cur = tok
                parts[cur] = []
            elif cur:
                parts[cur].append(tok)
        if parts.get("NAME"):
            self.tests[parts["NAME"][0]] = {"prefix": parts.get("PREFIX", []), "command": parts.get("COMMAND", [])}


def gust_run(test):
    """(deck, overrides, patterns, expects_abort, problem) of a test that runs gust, or None."""
    cmd, prefix = test["command"], test["prefix"]
    if GUST not in cmd:
        return None
    at = cmd.index(GUST)
    patterns, aborts, problem = [], False, ""
    if prefix and any(p.endswith("expect_abort.py") for p in prefix):
        k = next(i for i, p in enumerate(prefix) if p.endswith("expect_abort.py")) + 1
        aborts = True
        while k < len(prefix) and prefix[k] != "--":
            if prefix[k] == "--success":
                aborts = False
            elif prefix[k] in ("--log", "--timeout"):
                k += 1
            else:
                patterns.append(prefix[k])
            k += 1
    elif prefix:
        script = os.path.basename(prefix[1] if len(prefix) > 1 else prefix[0])
        if script not in TRANSPARENT:
            problem = "run by %s, which adds keys" % script
    if at != 0:
        problem = "gust is run by " + os.path.basename(cmd[1] if len(cmd) > 1 else cmd[0])
    rest = cmd[at + 1:]
    if not rest or "=" in rest[0]:
        return None, [], patterns, aborts, problem or "no deck on the command line"
    deck, overrides = rest[0], rest[1:]
    if any("${" in x or "$<" in x for x in [deck] + overrides):
        problem = problem or "an unresolved variable"
    elif not os.path.isfile(deck):
        problem = problem or "the deck is written at run time"
    elif any("=" not in x for x in overrides):
        problem = problem or "an argument that is not key=value"
    return deck, overrides, patterns, aborts, problem


def named_keys(pattern):
    """The GUST keys a regex of expect_abort.py names."""
    text = pattern.replace("\\.", ".")
    return [m for m in re.findall(r"[a-z]+\.[A-Za-z0-9_]+", text) if m.split(".")[0] in OWNED]


# ---------------------------------------------------------------------------
# The tests
# ---------------------------------------------------------------------------
# Decks that GUST refuses on a rule the checker reports under another key of the same rule, or that only the run
# can see (a restart checkpoint's contents): the checker reports an ERROR, but not on the key the message names.
OTHER_KEY = {
    "abort_force_file_restart_header",        # io.restart: the checkpoint is written by the test's first run
    "abort_unsteady_restart_average_start",   # the same
    "abort_unsteady_restart_past_stop",       # the same
    "abort_buildings_lusgs",                  # reported on steady.sweep
    "abort_abl_capped_no_sink",               # reported on bc.rayleigh_fields / bc.theta_profile
    "abort_stage_no_guess",                   # reported on init.first_guess
    "abort_sequencing_overlap",               # reported on steady.line_overlap
}
MIN_STARTING = 90     # decks of the CTests that must start and that the sweep can resolve (98 in November 2026)
MIN_ABORTING = 95     # decks that must stop and whose message names a key (105 in November 2026)


def ctest_decks():
    tests = CMakeTests(os.path.join(TESTS, "CMakeLists.txt"))
    for name, test in tests.tests.items():
        run = gust_run(test)
        if run is None or run[4]:
            continue
        deck, overrides, patterns, aborts, _ = run
        yield name, deck, overrides, patterns, aborts


def errors(deck, overrides):
    report = gic.check_deck(deck, overrides=tuple(overrides), options=Options(cwd=os.path.dirname(deck)))
    return [f for f in report.findings if f.severity == Severity.ERROR]


class SelfTestTable(unittest.TestCase):
    def test_every_case(self):
        for case in selftest.CASES:
            with self.subTest(case=case[0]):
                problems, _ = selftest.run_case(case)
                self.assertEqual(problems, [])

    def test_case_names_unique(self):
        names = [case[0] for case in selftest.CASES]
        self.assertEqual(len(names), len(set(names)))

    def test_every_checker_has_a_tripping_case(self):
        from gust_input_checker.checker import registry
        missing = set(c.name for c in registry()) - selftest.checkers_covered()
        self.assertEqual(missing, set())


class ValueTable(unittest.TestCase):
    def test_every_reference_key_has_a_row(self):
        ref = load_reference(None)
        missing = sorted(k for k in ref.keys if k.split(".")[0] in OWNED and k not in values.SPEC)
        self.assertEqual(missing, [])

    def test_no_row_for_a_key_gust_does_not_read(self):
        ref = load_reference(None)
        self.assertEqual(sorted(k for k in values.SPEC if not ref.has(k)), [])

    def test_tke_floor_range(self):
        # Closure.cpp: 0 or at least machine epsilon.
        row = values.SPEC["turb.tke_floor"]
        self.assertEqual(row.range, "0|>=2.220446049250313e-16")


class CTestDecks(unittest.TestCase):
    """The decks of Tests/CMakeLists.txt, with the keys each test adds on the command line."""

    @classmethod
    def setUpClass(cls):
        cls.decks = list(ctest_decks())

    def test_decks_that_start_have_no_error(self):
        count = 0
        for name, deck, overrides, _, aborts in self.decks:
            if aborts:
                continue
            count += 1
            with self.subTest(test=name):
                # A file the test writes before GUST runs (or a checkpoint of an earlier test) is not there yet.
                found = [(f.key, f.message) for f in errors(deck, overrides) if f.checker != "files"]
                self.assertEqual(found, [])
        self.assertGreaterEqual(count, MIN_STARTING)

    def test_decks_that_stop_have_an_error_on_the_named_key(self):
        count = 0
        for name, deck, overrides, patterns, aborts in self.decks:
            keys = set(k for p in patterns if not p.startswith("!") for k in named_keys(p))
            if not aborts or not keys:
                continue
            count += 1
            with self.subTest(test=name):
                found = errors(deck, overrides)
                self.assertTrue(found, "no ERROR for a deck GUST refuses")
                on_key = any(f.key in keys for f in found)
                self.assertEqual(on_key, name not in OTHER_KEY,
                                 "ERRORs on %s; GUST names %s" % (sorted(set(f.key for f in found)), sorted(keys)))
        self.assertGreaterEqual(count, MIN_ABORTING)


class CommandLine(unittest.TestCase):
    def run_tool(self, *args):
        import subprocess
        import sys
        return subprocess.run([sys.executable, os.path.join(HERE, "check_gust_inputs.py")] + list(args),
                              capture_output=True, text=True)

    def test_exit_status(self):
        import tempfile
        with tempfile.TemporaryDirectory() as folder:
            deck = os.path.join(folder, "inputs")
            with open(deck, "w") as fh:
                fh.write(selftest.SOLVE)
            self.assertEqual(self.run_tool(deck, "--min-severity", "ERROR").returncode, 0)
            bad = self.run_tool(deck, "steady.cfl_max=-1")
            self.assertEqual(bad.returncode, 1)
            self.assertIn("steady.cfl_max", bad.stdout)
            made = self.run_tool(deck, "init.first_guess=column", "bc.inflow_profile=file",
                                 "bc.inflow_file=column/t.dat")
            self.assertIn("bc.inflow_file", made.stdout)
            later = self.run_tool(deck, "bc.inflow_profile=file", "bc.inflow_file=column/t.dat",
                                  "--made-by-run", "column/t.dat", "--min-severity", "ERROR")
            self.assertNotIn("[ERROR] bc.inflow_file", later.stdout)


if __name__ == "__main__":
    unittest.main()
