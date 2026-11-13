#ifndef GUST_UNSTEADY_OUTPUT_H_
#define GUST_UNSTEADY_OUTPUT_H_

// The outputs of an unsteady run (steady.stage = unsteady, TimeDriver.H):
// running time averages and variances, probe time series, and the record a
// checkpoint needs so that a restart continues the run.
//
// Averages (RunningAverages). After every physical step t^n -> t^{n+1} the
// cell-centred fields x = (u_c, v_c, w_c, theta, k) at t^{n+1} are added with
// the weight
//
//   w = t^{n+1} - max(t^n, unsteady.average_start)      (no sample if w <= 0),
//
// the length of the step's part after the averaging start. With BDF1, q^{n+1}
// stands for the whole step, so this is the rectangle rule of the time
// integral. The weighted update of West (1979), with W the sum of the weights,
//
//   W <- W + w,   d = x - mean,   mean <- mean + (w/W) d,   M2 <- M2 + w d (x - mean),
//
// keeps the mean and the variance M2/W without the cancellation of a sum of
// squares. The variance is that of the resolved unsteadiness; the modelled
// part is k. u_c, v_c and w_c are the face values averaged to the cell centre
// (as in the plotfile).
//
// Probes (Probes). unsteady.probes_file lists points, one per line,
// "x y h" with h the height above the ground [m]; '#' starts a comment. Each
// probe takes the cell whose column holds (x, y) and whose centre height
// above the column's ground (the mean of its four ground nodes) is nearest h
// (the lower cell on a tie). No interpolation: the cell's values. After every
// step one row "t, u_c, v_c, w_c, theta, k for each probe" goes to
// <io.plot_file>_probes.csv, written by the I/O rank; a header of '#' lines
// names every probe's point and cell.
//
// The checkpoint record (write_unsteady_record / read_unsteady_record):
// <checkpoint>/Unsteady (text, hex floats: the time, the steps, the fixed
// step or 0 with unsteady.cfl, unsteady.average_start and W) and the
// averages' mean and M2 as VisMF MultiFabs (Level_0/unsteady_mean,
// Level_0/unsteady_m2). q^n is not written: at the end of a step it is the
// state itself, and the next step stores it (TimeTerm::store).

#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! The fields averaged and probed, in order: u, v, w, theta, tke.
inline constexpr int kUnsteadyFields = 5;
const std::vector<std::string>& unsteady_field_names ();

class RunningAverages
{
public:
    //! On (ba, dm), the state's cell-centred layout; nothing averaged yet.
    RunningAverages (amrex::Real average_start, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! Adds q (at t1, the step [t0, t1]) with the weight t1 - max(t0, average_start), when > 0.
    //! Reads the valid faces and cells of q only.
    void add (const State& q, double t0, double t1);

    [[nodiscard]] double weight () const { return m_weight; }
    [[nodiscard]] amrex::Real average_start () const { return m_start; }
    [[nodiscard]] const amrex::MultiFab& mean () const { return m_mean; }
    [[nodiscard]] const amrex::MultiFab& m2 () const { return m_m2; }

    //! The plotfile's averaged fields: u_mean ... tke_mean, then u_var ... tke_var (M2/W;
    //! 0 everywhere while W = 0), and their names.
    void plot_fields (amrex::MultiFab& out, std::vector<std::string>& names) const;

    //! The checkpoint's part: <dir>/Level_0/unsteady_mean and unsteady_m2.
    void write_checkpoint (const std::string& dir) const;
    //! Restores W, the mean and M2 (a checkpoint of another box layout or rank count is fine).
    void read_checkpoint (const std::string& dir, double weight);

private:
    amrex::Real m_start;
    double m_weight{0};
    amrex::MultiFab m_mean, m_m2;
};

//! One probe: the point as given and the cell it takes.
struct Probe
{
    amrex::Real x{0}, y{0}, h{0};   //!< as given [m]
    amrex::IntVect cell;            //!< the cell's index
    amrex::Real h_cell{0};          //!< the cell centre's height above the ground [m]
};

//! Reads unsteady.probes_file's text: "x y h" per line, '#' comments, blank lines skipped.
//! Errors (naming unsteady.probes_file and the line) go to err; no probe is an error.
std::vector<Probe> parse_probes (const std::string& text, const std::string& file, std::vector<std::string>& err);

class Probes
{
public:
    //! Locates every probe on the mesh (collective): the cell centre height is the mean of the
    //! cell's eight node heights, as the plotfile's z. Aborts, naming probes_file
    //! (unsteady.probes_file), on a point outside the domain in x or y, or a height above the
    //! column's top.
    Probes (std::vector<Probe> probes, const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom,
            const std::string& probes_file);

    [[nodiscard]] const std::vector<Probe>& probes () const { return m_probes; }

    //! The values at the probes: kUnsteadyFields per probe, the same on every rank (collective).
    std::vector<amrex::Real> sample (const State& q) const;

    //! The '#' header lines and the column line of the output file.
    [[nodiscard]] std::string header () const;

    //! Starts the output file out (the I/O rank writes; collective): a new file with the header;
    //! or, with keep_until >= 0 (a restart from an unsteady checkpoint at that time), the
    //! existing file kept up to the rows with t <= keep_until and appended to. Returns true when
    //! it kept an existing file, false when it started a new one (the caller writes the first
    //! row); prints a notice when a restart finds no file or a new run replaces one. Aborts,
    //! naming unsteady.probes_file, when the existing file's header is not this one.
    bool open (const std::string& out, const std::string& probes_file, double keep_until) const;
    //! Appends one row (I/O rank writes; collective, as sample()).
    void write_row (const std::string& out, double t, const State& q) const;

private:
    bool open_io (const std::string& out, const std::string& probes_file, double keep_until) const;
    std::vector<Probe> m_probes;
};

//! What a checkpoint of an unsteady run records (besides the averages' fields).
struct UnsteadyRecord
{
    double time{0};          //!< t^n [s]
    long steps{0};           //!< physical steps run to t^n
    double dt{0};            //!< the fixed step unsteady.dt [s]; 0 with unsteady.cfl
    double average_start{0}; //!< unsteady.average_start [s]
    double weight{0};        //!< W, the averages' summed weight [s]
    double anchor_time{0};   //!< the fixed step's time anchor t_* [s] (TimeDriver.H)
    long anchor_steps{0};    //!< its total step n_a
};

//! <dir>/Unsteady, text with hex floats (I/O rank).
void write_unsteady_record (const std::string& dir, const UnsteadyRecord& r);
//! true and r filled when <dir>/Unsteady exists (a checkpoint of an unsteady run); false for a
//! steady checkpoint. Aborts, naming io.restart, on a file that does not parse. Collective.
bool read_unsteady_record (const std::string& dir, UnsteadyRecord& r);

//! The probe output file of a run: <io.plot_file>_probes.csv.
std::string probes_output_file (const std::string& plot_file);

} // namespace gust

#endif
