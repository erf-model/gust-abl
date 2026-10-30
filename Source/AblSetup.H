#ifndef GUST_ABLSETUP_H_
#define GUST_ABLSETUP_H_

// Named atmospheric boundary layer set-ups: abl.drive (what pushes the flow)
// and abl.theta (how theta is treated). Neither key adds physics. Each name
// gives a few other keys a default when the deck leaves them out, and then
// checks that the deck's final inputs are that set-up; a key that contradicts
// the name stops the run, naming both. Without the keys nothing changes.
//
//   abl.drive            defaults set                     needs
//   inflow_only          -                                inflow face, no forcing
//   couette              bc.inflow_tke_length = unbounded log-law inflow, moving no-slip lid,
//                        turb.max_geom_lscale = 1e6       no forcing, abl.theta = neutral
//   pressure_driven      -                                pressure gradient, slip top, no Coriolis
//   geostrophic_column   forcing.coriolis = 1,            a geostrophic or column target wind,
//                        forcing.geostrophic = 1          column first guess, inflow faces take a table
//   geostrophic_inflow   as geostrophic_column            a geostrophic wind, inflow face
//
//   abl.theta            defaults set                     needs
//   neutral              -                                uniform theta, neutral surface, no leaf heat
//   capped               bc.rayleigh_fields = w theta     rising sounding, top gradient, Rayleigh layer on w and theta
//   frozen               simplec.theta = frozen           sounding or inflow table, neutral surface
//   surface_flux         bc.surface_stability = heat_flux surface heat flux, and a heat sink (outflow face,
//                                                         Rayleigh layer on theta or a marching column)
//
// With a file inflow, theta means the table's. pressure_driven also needs: no
// Rayleigh layer on u, v over uniform theta, one over a stratified theta.

#include "BoundaryConditions.H"
#include "Closure.H"
#include "Forcing.H"
#include "Initializer.H"
#include "Simplec.H"
#include "SurfaceLayer.H"

#include <string>
#include <vector>

namespace gust {

struct AblSetupOptions
{
    enum class Drive { unset, inflow_only, couette, pressure_driven, geostrophic_column, geostrophic_inflow };
    enum class Theta { unset, neutral, capped, frozen, surface_flux };

    Drive drive{Drive::unset};
    Theta theta{Theta::unset};
    //! Bad names and repeated keys, reported with the run's other input errors.
    std::vector<std::string> errors;

    [[nodiscard]] bool active () const { return drive != Drive::unset || theta != Theta::unset; }
};

//! The deck's inputs that a set-up is checked against. Value copies of what
//! the application read, so that the check needs no ParmParse.
struct AblSetupFacts
{
    BoundaryInputs bc;
    ForcingOptions forcing;
    InitOptions init;
    SurfaceLayerOptions surface;
    SimplecOptions simplec;
    //! theta of the bc.inflow_file table, with bc.inflow_profile = file.
    amrex::Vector<amrex::Real> inflow_table_theta;
    amrex::Real max_geom_lscale{30};
    bool surface_layer{false};  //!< bc.zlo = surface_layer
    bool restart{false};        //!< io.restart is given (the first guess is the checkpoint)
    bool solve{false};          //!< steady.stage is not none
    bool canopy_heat{false};    //!< canopy.heat is on
};

//! The drive's and the theta treatment's input names (for messages and logs).
const char* abl_drive_name (AblSetupOptions::Drive d);
const char* abl_theta_name (AblSetupOptions::Theta t);

//! Read abl.drive and abl.theta; a bad name, a repeated key or a deck value
//! that contradicts a name's default goes to errors, which the application
//! reports before any other reader.
AblSetupOptions read_abl_setup_options ();
//! Give the set-up's keys their defaults, in ParmParse, where the deck leaves
//! them out. Called before any reader; does nothing without a stage (the keys
//! are then refused by check_abl_setup) or after a name error.
void apply_abl_setup_defaults (const AblSetupOptions& opt);
//! Every way the deck's final inputs differ from the set-up, each message
//! naming abl.drive or abl.theta and the input.
std::vector<std::string> check_abl_setup (const AblSetupOptions& opt, const AblSetupFacts& facts);

} // namespace gust
#endif
