"""Facts derived once from a deck, shared by every check.

`RunContext` answers "what does GUST do with this deck": which stage runs,
which readers run (and so which keys are read at all), the faces, the vertical
levels and their top, the defaults a named set-up (abl.drive, abl.theta) adds.
The accessors never raise: a value that cannot be read gives None, and the
`values` check reports it.

Sources: Source/App.cpp (run: which readers run in which mode),
Source/AblSetup.cpp (apply_abl_setup_defaults), Source/ZLevels.cpp
(make_zlevels) and Tests/inputs_reference (the defaults).
"""
import fnmatch
import math
import os

from .deck import Unparsed, to_int, to_real
from .reference import load_reference

FACES = ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi")
LATERAL = FACES[:4]
OWNED = ("mesh", "terrain", "init", "turb", "bc", "steady", "solver", "io", "forcing", "convection", "simplec",
         "canopy", "buildings", "abl", "unsteady", "nest")
AMREX = ("amrex", "fab", "fabarray", "vismf", "tiny_profiler", "DistributionMapping", "device", "hypre")
SURFACE_KEYS = ("surface_z0", "surface_z0_file", "surface_zref", "surface_stability", "surface_heat_flux",
                "surface_temperature", "surface_temperature_file", "surface_kb_inverse",
                "surface_wstar", "surface_pblh", "surface_pblh_file",
                "surface_wall_k", "surface_denominator_floor", "surface_stable_functions", "surface_qnse_am",
                "surface_qnse_bm", "surface_qnse_ah", "surface_qnse_bh")
NEWTON_KEYS = ("newton_gmres_restart", "newton_gmres_max_iters", "newton_fd_eps_scale", "newton_krylov_rtol",
               "newton_mem_budget_gb", "newton_armijo_c", "newton_backtrack", "newton_lambda_min", "newton_cfl_cut",
               "newton_max_retries", "newton_dtau_inf_at_cap", "newton_ptc_growth")
SEQUENCING_KEYS = ("sequencing", "sequencing_drop_orders", "sequencing_monitor_tol", "sequencing_hold",
                   "sequencing_cfl_start")
RAYLEIGH_FIELDS = ("u", "v", "w", "theta", "tke")


class Options:
    """What the user says about the run that the deck does not hold."""

    def __init__(self, ranks=None, gpu=False, reference=None, cwd=None, made_by_run=()):
        self.ranks = ranks              # MPI ranks of the planned run (None: not said)
        self.gpu = gpu                  # the run uses a GPU build
        self.reference = reference      # path of Tests/inputs_reference (None: this repository's)
        self.cwd = cwd                  # the directory the run starts in (None: the current one)
        # Paths or glob patterns of files and directories that an earlier step of the run writes: a key that names
        # one of them is not reported missing.
        self.made_by_run = (made_by_run,) if isinstance(made_by_run, str) else tuple(made_by_run or ())


def geometric_ratio(dz0, top, n):
    """The ratio r with dz0 (1 + r + ... + r^(n-1)) = top, or None (ZLevels.cpp: stretch to a fixed top)."""
    if not (dz0 > 0 and top > 0 and n >= 1) or dz0 * n > top * (1 + 1e-12):
        return None
    if n == 1 or abs(dz0 * n - top) <= 1e-12 * top:
        return 1.0
    lo, hi = 1.0, 2.0

    def total(r):
        return dz0 * n if r == 1.0 else dz0 * (r ** n - 1.0) / (r - 1.0)

    while total(hi) < top and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < top:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def gust_newton_converges(dz0, top, n, tol=1.0e-8):
    """Whether GUST finds the ratio itself (mesh.dz0 without mesh.stretch_ratio).

    A copy of ERF's geometric_stretch (extern_erf/ERF_GridUtils.H), which GUST
    calls in double precision with tol = 1e-8 (ZLevels.cpp): Newton from
    r = 1.03, at most 50 steps, r kept at or above 1 + tol. It fails for steep
    stretches (r^n beyond about 1e3) where the bisection of geometric_ratio
    still finds the ratio, and just below dz0 = top / n.
    """
    if not (dz0 > 0 and top > 0 and n >= 1):
        return False
    if n == 1 or dz0 >= top / n:
        return True

    def residual(r):
        return dz0 * ((r ** n - 1.0) / (r - 1.0)) - top

    r = 1.03
    try:
        f = residual(r)
        steps = 0
        while abs(f) > tol and steps < 50:
            dfdr = dz0 * (n * r ** (n - 1.0) * (r - 1.0) - r ** n + 1.0) / (r - 1.0) ** 2.0
            r -= f / dfdr
            r = max(r, 1.0 + tol)
            f = residual(r)
            steps += 1
    except (OverflowError, ZeroDivisionError):
        return False
    return abs(f) <= tol


class RunContext:
    """Derived facts of one deck."""

    def __init__(self, deck, options=None):
        self.deck = deck
        self.options = options or Options()
        self.problems = []          # (key, message): reported as ERROR by run_all_checks
        self.ref = load_reference(self.options.reference)
        self.cwd = self.options.cwd or os.getcwd()
        self.deck_dir = os.path.dirname(os.path.abspath(deck.path)) if deck.path else self.cwd
        self.abl_added = {}
        self._stage_and_abl()
        self._mesh()
        self._boundaries()

    # -- access -----------------------------------------------------------
    def given(self, key):
        """Whether GUST's readers see the key: in the deck, or added by a named set-up (abl.*)."""
        return self.deck.has(key) or key in self.abl_added

    def tokens(self, key):
        if self.deck.has(key):
            return self.deck.values(key)
        if key in self.abl_added:
            return list(self.abl_added[key])
        if key == "steady.sweep":
            return ["zjacobi" if self.options.gpu else "lusgs"]
        if key == "steady.cfl_max":
            return ["10" if self.s("steady.sweep") == "zjacobi" else "8"]
        if key == "steady.newton_fd_eps_scale":
            return [repr(math.sqrt(2.220446049250313e-16))]
        return self.ref.default(key)

    def s(self, key):
        """The first value as a string ('' when the key has no value)."""
        toks = self.tokens(key)
        return toks[0] if toks else ""

    def rs(self, key):
        """The values as numbers, or None when one of them cannot be read."""
        out = []
        for t in self.tokens(key):
            try:
                out.append(to_real(t))
            except (ValueError, Unparsed):
                return None
        return out

    def r(self, key):
        vals = self.rs(key)
        return vals[0] if vals else None

    def i(self, key):
        toks = self.tokens(key)
        if not toks:
            return None
        try:
            return to_int(toks[0])
        except (ValueError, Unparsed):
            return None

    def ints(self, key):
        out = []
        for t in self.tokens(key):
            try:
                out.append(to_int(t))
            except (ValueError, Unparsed):
                return None
        return out

    def flag(self, key):
        return self.i(key) == 1

    def positive(self, key):
        x = self.r(key)
        return x is not None and x > 0

    # -- the stage and the named set-ups ---------------------------------
    def _stage_and_abl(self):
        self.stage = self.s("steady.stage")
        self.solve = self.stage in ("simplec", "simplec_newton", "unsteady")
        self.newton = self.stage == "simplec_newton"
        self.unsteady = self.stage == "unsteady"
        self.first_guess = self.s("init.first_guess")
        self.initialize = self.first_guess in ("profile", "column", "nest")
        self.reads_bc = self.initialize or self.solve      # bc.* and solver.* are read (Source/App.cpp)
        self.restart = self.s("io.restart") if self.deck.has("io.restart") else ""
        self.abl_drive = self.s("abl.drive") if self.deck.has("abl.drive") else ""
        self.abl_theta = self.s("abl.theta") if self.deck.has("abl.theta") else ""
        drives = ("inflow_only", "couette", "pressure_driven", "geostrophic_column", "geostrophic_inflow")
        thetas = ("neutral", "capped", "frozen", "surface_flux")
        self.abl_valid = ((not self.abl_drive or self.abl_drive in drives)
                          and (not self.abl_theta or self.abl_theta in thetas)
                          and self.deck.count("abl.drive") <= 1 and self.deck.count("abl.theta") <= 1)
        # AblSetup.cpp, read_abl_setup_options: a deck value the name refuses (compared as text, as the source does).
        self.abl_refused = []
        rules = []
        if self.abl_valid and self.abl_drive in ("geostrophic_column", "geostrophic_inflow"):
            rules += [("forcing.coriolis", "0", "abl.drive = " + self.abl_drive),
                      ("forcing.geostrophic", "0", "abl.drive = " + self.abl_drive)]
        if self.abl_valid and self.abl_theta == "frozen":
            rules.append(("simplec.theta", "solve", "abl.theta = frozen"))
        if self.abl_valid and self.abl_theta == "surface_flux":
            rules.append(("bc.surface_stability", "neutral", "abl.theta = surface_flux"))
        for key, bad, name in rules:
            if self.deck.has(key) and self.deck.values(key)[:1] == [bad]:
                self.abl_refused.append((key, "%s refuses %s = %s" % (name, key, bad)))
        # AblSetup.cpp, apply_abl_setup_defaults: only with a stage and no abl error; a key of the deck is kept.
        if self.stage != "none" and self.abl_valid and not self.abl_refused:
            add = {}
            if self.abl_drive in ("geostrophic_column", "geostrophic_inflow"):
                add["forcing.coriolis"] = ["1"]
                add["forcing.geostrophic"] = ["1"]
            if self.abl_drive == "couette":
                add["bc.inflow_tke_length"] = ["unbounded"]
                add["turb.max_geom_lscale"] = ["1e6"]
            if self.abl_theta == "capped":
                add["bc.rayleigh_fields"] = ["w", "theta"]
            if self.abl_theta == "frozen":
                add["simplec.theta"] = ["frozen"]
            if self.abl_theta == "surface_flux":
                add["bc.surface_stability"] = ["heat_flux"]
            self.abl_added = {k: v for k, v in add.items() if not self.deck.has(k)}

    # -- the mesh ---------------------------------------------------------
    def _mesh(self):
        n = self.ints("mesh.n_cell")
        self.n_cell = n if n and len(n) == 3 and all(v >= 1 for v in n) else None
        lo, hi = self.rs("mesh.prob_lo"), self.rs("mesh.prob_hi")
        self.prob_lo = lo if lo and len(lo) == 3 else None
        self.prob_hi = hi if hi and len(hi) == 3 else None
        self.extent = None
        if self.prob_lo and self.prob_hi and all(h > l for l, h in zip(self.prob_lo, self.prob_hi)):
            self.extent = [h - l for l, h in zip(self.prob_lo, self.prob_hi)]
        per = self.ints("mesh.is_periodic")
        self.periodic = per if per and len(per) == 3 else None
        self.dx = self.dy = None
        if self.n_cell and self.extent:
            self.dx = self.extent[0] / self.n_cell[0]
            self.dy = self.extent[1] / self.n_cell[1]
        mgs = self.ints("mesh.max_grid_size")
        if mgs and len(mgs) in (1, 2) and all(v >= 1 for v in mgs):
            self.max_grid = (mgs[0], mgs[-1])
        else:
            self.max_grid = None
        lev = self.i("mesh.hierarchy_levels")
        self.levels = lev if lev is not None and 1 <= lev <= 10 else 1
        self.terrain_type = self.s("terrain.type")
        self.has_terrain = self.terrain_type in ("agnesi", "xyz", "erf_grid")
        self._zlevels()

    def _zlevels(self):
        """The node heights of a flat column (ZLevels.cpp), where the deck decides them."""
        self.zlevels = None
        self.zlevels_source = ""
        self.stretch = None            # the ratio between successive layers
        if not (self.n_cell and self.extent):
            return
        nz, top = self.n_cell[2], self.extent[2]
        dz0, ratio = self.r("mesh.dz0"), self.r("mesh.stretch_ratio")
        if self.deck.has("mesh.zlevels_file"):
            self.zlevels_source = "file"
            found = self.find_file(self.s("mesh.zlevels_file"))[0]
            if found:
                vals = read_numbers(found)
                if vals and len(vals) == nz + 1 and all(b > a for a, b in zip(vals, vals[1:])):
                    self.zlevels = vals
            return
        if dz0 is None or not dz0 > 0:
            self.zlevels_source = "uniform"
            self.zlevels = [k * top / nz for k in range(nz + 1)]
            self.stretch = 1.0
            return
        if self.deck.has("mesh.stretch_ratio"):
            self.zlevels_source = "ratio"
            if ratio is None or not (ratio >= 1 and math.isfinite(ratio)):
                return
            r = ratio
        else:
            self.zlevels_source = "fixed_top"
            r = geometric_ratio(dz0, top, nz)
            if r is None:
                return
        z, dz = [0.0], dz0
        for _ in range(nz):
            z.append(z[-1] + dz)
            dz *= r
        if not math.isfinite(z[-1]):
            return
        if self.zlevels_source == "fixed_top":
            z[-1] = top
        self.zlevels = z
        self.stretch = r

    @property
    def ztop(self):
        return self.zlevels[-1] if self.zlevels else None

    @property
    def dz_first(self):
        return self.zlevels[1] - self.zlevels[0] if self.zlevels else None

    @property
    def dz_last(self):
        return self.zlevels[-1] - self.zlevels[-2] if self.zlevels else None

    # -- faces and what the stage reads ----------------------------------
    def _boundaries(self):
        self.face = {f: self.s("bc." + f) for f in FACES}
        self.surface_layer = self.reads_bc and self.face["zlo"] == "surface_layer"
        self.inflow_faces = [f for f in FACES if self.face[f] == "inflow"]
        self.outflow_faces = [f for f in FACES if self.face[f] == "outflow"]
        self.nest_faces = [f for f in LATERAL if self.face[f] == "nest"]
        self.inflow_profile = self.s("bc.inflow_profile")
        # App.cpp: the profile is needed by the first guess (not a nest's) and by a restart.
        self.profile_needed = (self.initialize and self.first_guess != "nest") or bool(self.restart)
        self.needs_profile = (bool(self.inflow_faces) or self.positive("bc.sponge_width")
                              or self.positive("bc.rayleigh_zdamp") or self.profile_needed)
        self.coriolis = self.solve and self.flag("forcing.coriolis")
        self.geostrophic = self.coriolis and self.flag("forcing.geostrophic")
        self.closure = self.s("turb.closure")
        self.realizability = self.r("turb.realizability") if self.solve else 0.0
        sig = self.r("convection.hybrid_sigma_max")
        self.hybrid = self.solve and sig is not None and sig > 0
        self.theta_frozen = self.solve and self.s("simplec.theta") == "frozen"
        self.sweep = self.s("steady.sweep")
        self.sequencing = self.solve and self.flag("steady.sequencing")
        self.canopy = self.s("canopy.type") in ("patches", "map")
        self.buildings = self.s("buildings.type") == "immersed_forcing"
        theta = self.rs("bc.theta_profile")
        self.stratified = bool(theta) and len(theta) > 1 and any(t != theta[0] for t in theta)

    @property
    def heat_flux_surface(self):
        """bc.surface_stability = heat_flux on a surface-layer ground."""
        return self.surface_layer and self.s("bc.surface_stability") == "heat_flux"

    @property
    def temperature_surface(self):
        """bc.surface_stability = surface_temperature on a surface-layer ground."""
        return self.surface_layer and self.s("bc.surface_stability") == "surface_temperature"

    @property
    def column_evolves(self):
        """Whether the column first guess solves theta (Source/ColumnGuess.cpp, column_theta_evolves)."""
        mode = self.s("init.column_stability")
        if mode == "mol":
            mol = self.r("init.column_mol")
            return mol is not None and mol != 0
        if mode == "heat_flux":
            return self.heat_flux_surface
        return mode == "surface_temperature"

    @property
    def aspect_ratio(self):
        """max(dx, dy) / first layer: the cell aspect ratio over flat ground (terrain only raises it), or None."""
        first = self.dz_first
        if not self.dx or not self.dy or not first or first <= 0:
            return None
        return max(self.dx, self.dy) / first

    def rayleigh_field(self, name):
        """Whether bc.rayleigh_fields holds `name` (default: all five)."""
        return name in self.tokens("bc.rayleigh_fields")

    def rayleigh_damps(self, name):
        return self.positive("bc.rayleigh_zdamp") and self.positive("bc.rayleigh_dampcoef") and self.rayleigh_field(name)

    @property
    def coriolis_parameter(self):
        """f = 4 pi sin(latitude) / rotation_period (Source/Forcing.cpp), or None."""
        lat, period = self.r("forcing.latitude"), self.r("forcing.rotation_period")
        if lat is None or period is None or not period > 0:
            return None
        return 4.0 * math.pi * math.sin(math.radians(lat)) / period

    # -- files ------------------------------------------------------------
    def find_file(self, name):
        """(path, where) of a file a key names: GUST opens it relative to the directory the run starts in.

        `where` is 'cwd' when it is found there, 'deck' when it is found only
        beside the deck, '' when it is not found.
        """
        if not name:
            return None, ""
        first = name if os.path.isabs(name) else os.path.join(self.cwd, name)
        if os.path.exists(first):
            return first, "cwd"
        beside = os.path.join(self.deck_dir, name)
        if not os.path.isabs(name) and os.path.exists(beside):
            return beside, "deck"
        return None, ""


    def made_by_run(self, name):
        """Whether `name` (a key's value) is one of Options.made_by_run: as written, or as the same path seen from
        the deck's directory or from the directory the run starts in."""
        def forms(text):
            out = {text, os.path.normpath(text)}
            if not os.path.isabs(text):
                out.add(os.path.normpath(os.path.join(self.deck_dir, text)))
                out.add(os.path.normpath(os.path.join(self.cwd, text)))
            return out

        mine = forms(name)
        for pattern in self.options.made_by_run:
            theirs = forms(os.fspath(pattern))
            if mine & theirs or any(fnmatch.fnmatchcase(a, b) for a in mine for b in theirs):
                return True
        return False


def read_numbers(path, limit=2000000):
    """Every number of a small text file ('#' comments), or None when a token is not a number."""
    out = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                for tok in line.split("#", 1)[0].split():
                    out.append(float(tok))
                    if len(out) > limit:
                        return None
    except (OSError, ValueError):
        return None
    return out
