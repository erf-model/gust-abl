"""Keys that GUST does not read.

Covers, with the abort that backs each (docs/inputs.rst, "The run stops
before doing any work, naming the input, when ..."):

- a key without a namespace, or in a namespace that is neither GUST's nor
  AMReX's (Source/Inputs.cpp, unknown_inputs: "inputs are namespace.key",
  "unknown namespace");
- a key that is not in Tests/inputs_reference (unknown_inputs: "not an input
  of namespace"), with the nearest known name;
- a key of GUST's namespaces given more than once, in the deck or again on the
  command line (every reader: "is given more than once"; docs/running.rst);
- bc.* and solver.* keys in a run that builds no initial field and has no
  stage: they are read only with init.first_guess = profile, column or nest,
  or with a stage (Source/App.cpp), so GUST reports them as unknown.
"""
import difflib

from ..checker import Checker, register
from ..context import AMREX, OWNED


@register
class KeysChecker(Checker):
    """Unknown, misspelt and repeated keys; keys of a namespace this run does not read."""

    name = "keys"

    def run(self):
        ref = self.ctx.ref
        known = [k for k in ref.keys if k.split(".", 1)[0] in OWNED]
        for key in self.deck.keys():
            if "." not in key:
                if key == "FILE":
                    pass                # read in place by the deck reader (a file it cannot open is a deck problem)
                else:
                    self.error(key, "GUST stops: unknown input '%s': inputs are namespace.key." % key,
                               self._suggest(key, known))
                continue
            ns = key.split(".", 1)[0]
            if ns in AMREX:
                continue
            if ns not in OWNED:
                if ns == "refavg":
                    self.error(key, "GUST stops: unknown namespace 'refavg' (the refavg.* keys are inputs of the "
                                    "gust_refavg tool, not of gust).", "Give the key to gust_refavg, in its own file.")
                else:
                    self.error(key, "GUST stops: unknown input '%s': unknown namespace '%s'." % (key, ns),
                               "GUST's namespaces are %s; AMReX's are %s." % (", ".join(OWNED), ", ".join(AMREX)))
                continue
            if not ref.has(key):
                self.error(key, "GUST stops: unknown input '%s': not an input of namespace '%s'." % (key, ns),
                           self._suggest(key, known))
                continue
            if ns in ("bc", "solver") and not self.ctx.reads_bc:
                self.error(key, "GUST stops: %s is not read in this run (the %s.* keys are read with "
                                "init.first_guess = profile, column or nest, or with steady.stage = simplec)."
                           % (key, ns), "Set init.first_guess or steady.stage, or remove the key.")
        for key, entries in self.deck.duplicates().items():
            if key.split(".", 1)[0] not in OWNED:
                continue
            places = ", ".join(e.where for e in entries)
            advice = "Keep one definition."
            if any(e.where == "command line" for e in entries) and any(e.where != "command line" for e in entries):
                advice = ("A key the file sets cannot be given again on the command line: edit the file, or remove "
                          "the key from it.")
            self.error(key, "GUST stops: %s is given more than once (%s)." % (key, places), advice)

    @staticmethod
    def _suggest(key, known):
        near = difflib.get_close_matches(key, known, n=2, cutoff=0.6)
        if not near and "." in key:
            ns, _, tail = key.partition(".")
            tails = {k.split(".", 1)[1]: k for k in known if k.startswith(ns + ".")}
            near = [tails[t] for t in difflib.get_close_matches(tail, list(tails), n=2, cutoff=0.6)]
        if not near:
            return "The inputs are listed in Tests/inputs_reference and docs/inputs.rst."
        return "Did you mean %s?" % " or ".join(near)
