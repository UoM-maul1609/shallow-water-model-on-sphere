# Spherical Shallow-Water Model

A parallel Fortran shallow-water model on a rotating sphere, designed primarily for idealised planetary-atmosphere experiments and Saturn polar-jet simulations.

The model solves the nonlinear shallow-water equations on a longitude–latitude grid, supports observed or idealised zonal-jet initial conditions, optional spherical metric terms, multiple Coriolis and subgrid closures, stochastic height perturbations, zonal-mean nudging, free-slip or radiative latitude boundaries, optional latitude sponge layers, MPI domain decomposition, and NetCDF output.

The current default configuration is aimed at a truncated northern Saturn domain rather than a pole-to-pole global grid.

---

## 1. Model equations

The prognostic variables are layer thickness `h` and horizontal velocity components `u` and `v`, where:

- `u` is eastward/zonal velocity;
- `v` is northward/meridional velocity;
- `h` is the active shallow-water layer thickness;
- `hs` is stationary lower-boundary/topographic height;
- `height = h + hs` is the total free-surface height.

In spherical longitude–latitude coordinates `(lambda, theta)`, the model represents the nonlinear rotating shallow-water system. In conservative form, the mass equation is of the form

```text
∂h/∂t + 1/(a cosθ) ∂(hu)/∂λ + 1/(a cosθ) ∂(hv cosθ)/∂θ = 0
```

where `a` is the planetary radius.

When spherical momentum metric terms are enabled, the momentum equations include the curvature terms associated with spherical coordinates, including contributions proportional to

```text
+ huv tanθ / a
- hu² tanθ / a
```

in the zonal and meridional momentum equations respectively.

The Coriolis parameter is

```text
f = 2 Ω sinθ
```

with

```text
Ω = 2π / rotation_period
```

and the gravity-wave speed is approximately

```text
c = sqrt(g h).
```

---

## 2. Numerical method

### 2.1 Horizontal discretisation

The main solver is a finite-difference / finite-volume-like Lax–Wendroff shallow-water scheme implemented in `lax_wendroff_ll` in `advection.f90`.

The longitude–latitude grid is regular in angular coordinates:

- longitude is periodic over `2π`;
- latitude spans the configured physical model domain;
- one halo cell is currently used around each MPI subdomain.

Grid metrics account for the shrinking zonal distance with latitude through factors involving `cos(theta)`.

### 2.2 Time integration

The dynamical timestep is controlled by

```fortran
nm1%dt
```

and the nominal number of steps is

```text
ntim = ceiling(runtime / dt).
```

The model does **not** currently provide automatic timestep adaptation.

The main shallow-water update uses the model's Lax–Wendroff predictor/corrector formulation. Optional Coriolis treatment is selected separately.

### 2.3 Coriolis schemes

```fortran
nm1%coriolis_scheme=0
```

uses the legacy/original Coriolis treatment.

```fortran
nm1%coriolis_scheme=1
```

uses a Crank–Nicolson-style time-centred Coriolis update.

For the spherical Saturn configuration, scheme 1 is generally the more appropriate option.

### 2.4 Spherical momentum metric terms

```fortran
nm1%momentum_metric_terms=0
```

disables the explicit spherical curvature terms and retains the legacy formulation.

```fortran
nm1%momentum_metric_terms=1
```

enables spherical momentum curvature terms and uses the corresponding gradient-wind treatment during initialisation.

For high-latitude planetary simulations, `momentum_metric_terms=1` is recommended.

---

## 3. Grid and domain

The global grid is controlled by

```fortran
nm1%ip
nm1%jp
```

where:

- `ip` is the number of longitude cells;
- `jp` is the number of latitude points.

Longitude covers a full periodic circle.

The requested latitude limits are

```fortran
nm1%slat
nm1%nlat
```

and are clipped by

```fortran
nm1%slat_thresh
nm1%nlat_thresh
```

using

```text
slat = max(slat, slat_thresh)
nlat = min(nlat, nlat_thresh).
```

The default Saturn example therefore requests `-90°` to `90°`, but actually runs over

```text
65°N to 86.5°N
```

because of the threshold values.

Latitude coordinates are constructed directly from the global latitude index. This avoids propagating latitude coordinates through MPI ranks and fixes the previous northern-halo `dthetan` error.

The latitude spacing is

```text
dtheta = (nlat - slat)/(jp - 1)
```

and the staggered latitude spacing currently uses

```text
dthetan = dtheta.
```

### Polar limitation

The longitude–latitude formulation contains terms involving `1/cos(theta)` and `tan(theta)`. The model should therefore **not** place ordinary grid-cell centres exactly at ±90°.

For polar studies, use a truncated domain such as the current `65°–86.5°N` configuration, or implement a dedicated polar treatment / alternative grid before attempting a true pole-to-pole simulation.

---

## 4. MPI parallelisation

The model uses MPI with a two-dimensional Cartesian decomposition.

```fortran
mp1%reorder=.true.
```

is retained.

MPI determines the processor layout using `MPI_Dims_create`, and the topology is created with:

- periodic longitude;
- non-periodic latitude.

All Cartesian-neighbour communication uses the Cartesian communicator, so rank reordering is supported correctly.

Halo exchanges explicitly pack and unpack edge data. This is important for east/west Fortran array columns, which are strided in memory and cannot safely be treated as contiguous MPI buffers without an MPI datatype or packing.

Example:

```bash
./run.sh 8 namelist.in
```

runs the model using eight MPI processes.

The process topology is printed at startup, for example:

```text
Cartesian topology: 4 2
```

---

## 5. Initial wind configurations

Two initial wind configurations are available.

### 5.1 Saturn wind profile

```fortran
nm1%initial_winds=1
```

reads an observed/reference Saturn zonal-wind profile from the NetCDF file specified by

```fortran
nm1%inputfile
```

The supplied example uses:

```fortran
nm1%inputfile='saturn_winds_vs_latitude_update.nc'
```

The input file is expected to contain:

```text
latitude(nlats)
wind(nlats)
```

The wind profile supports three transformations.

#### Wind amplitude

```fortran
nm1%wind_factor
```

multiplies the reference wind profile.

#### Latitude shift

```fortran
nm1%wind_shift
```

shifts profile features in latitude, in degrees. Positive values shift features northward.

#### Wind reduction

```fortran
nm1%wind_reduce
```

subtracts a constant velocity in m/s from the transformed wind profile.

Conceptually:

```text
Unew(theta) = wind_factor * Uref(theta - wind_shift) - wind_reduce.
```

A rotation-frame correction is then applied when the selected model rotation period differs from the reference System III period used by the wind data.

### 5.2 Idealised Gaussian jet

```fortran
nm1%initial_winds=2
```

creates an idealised Gaussian zonal jet using

```fortran
nm1%u_jet
nm1%theta_jet
nm1%h_jet
```

where:

- `u_jet` is peak jet speed in m/s;
- `theta_jet` is jet-centre latitude in degrees;
- `h_jet` is the Gaussian width in degrees.

---

## 6. Balanced initialisation

```fortran
nm1%initially_geostrophic=.true.
```

requests a wind field consistent with the initial height field after perturbations have been added.

The exact treatment depends on `momentum_metric_terms`.

With

```fortran
nm1%momentum_metric_terms=0
```

the legacy geostrophic treatment is used.

With

```fortran
nm1%momentum_metric_terms=1
```

the initialisation includes spherical curvature and uses a nonlinear/linearised gradient-wind treatment appropriate to the spherical momentum equations.

This is particularly important for high-latitude Saturn simulations, where curvature is not negligible.

---

## 7. Initial height perturbations

Random free-surface perturbations can seed jet instability:

```fortran
nm1%add_random_height_noise=.true.
```

Two schemes are available.

### 7.1 Legacy grid-cell noise

```fortran
nm1%height_noise_scheme=0
```

adds independent Gaussian perturbations to individual grid cells in the active jet latitude band.

This mode is resolution dependent and is retained mainly for reproducibility of older simulations.

### 7.2 Correlated physical-scale noise

```fortran
nm1%height_noise_scheme=1
```

generates a global random field and applies Gaussian smoothing in latitude and longitude so that the perturbation correlation scale is approximately resolution independent in physical distance.

Parameters are:

```fortran
nm1%height_noise_amplitude
nm1%height_noise_corr_length
```

where:

- `height_noise_amplitude` is the requested RMS height perturbation in metres;
- `height_noise_corr_length` is the Gaussian correlation length in metres.

The resulting perturbation is de-meaned and normalised over the active jet band before being applied.

For Saturn winds the perturbation band is currently approximately `75°–80°N`. For the idealised jet it is based on `theta_jet ± 3 h_jet`.

Because every MPI rank generates the same global random field using the same seed, the perturbation is independent of MPI decomposition.

---

## 8. Latitude boundary conditions

The latitude boundary treatment is selected with

```fortran
nm1%lat_boundary_scheme
```

### 8.1 Scheme 0 — legacy frozen exterior ghosts

```fortran
nm1%lat_boundary_scheme=0
```

retains the legacy fixed/frozen physical ghost-cell behaviour.

This is **not** a strict solid wall. Disturbance fluxes can interact with the fixed exterior state.

It is retained for backward comparison with older runs.

### 8.2 Scheme 1 — reflecting free-slip wall

```fortran
nm1%lat_boundary_scheme=1
```

implements a reflecting free-slip latitude wall.

The normal shallow-water mass flux is constrained to zero at the physical wall. Ghost values are constructed relative to the initialized balanced reference state so the background jet continuation is preserved as well as possible.

This option is useful for controlled closed-domain experiments, but it is generally **not physically appropriate for an artificial latitude cut through a planetary atmosphere**, because a real atmosphere has no rigid wall at, for example, 86.5°N.

Reflecting walls can also trap wave energy inside the domain.

### 8.3 Scheme 2 — balanced radiative/open boundary

```fortran
nm1%lat_boundary_scheme=2
```

is the recommended option for truncated planetary domains.

It operates on perturbations relative to the initialized balanced state.

For the linear normal shallow-water modes, the characteristic variables are approximately

```text
R± = v' ± g eta'/c
```

with

```text
c = sqrt(g href),
eta' = eta - eta_ref,
v' = v - v_ref.
```

At each physical latitude edge:

- the outgoing characteristic is taken from the adjacent interior row;
- the incoming perturbation characteristic is set to zero;
- the tangential wind perturbation uses zero normal gradient.

At the north edge, `R+` is treated as outgoing. At the south edge, `R-` is outgoing.

Unlike the free-slip wall, this condition does **not** impose zero meridional mass flux. It therefore allows disturbances to leave the truncated model domain.

The boundary is linearised around the initialized reference state. It is an approximate radiation condition rather than an exact non-reflecting boundary for all nonlinear, oblique, or Rossby-wave modes.

---

## 9. Latitude sponge layers

Optional sponge layers can be used with either:

```text
lat_boundary_scheme = 1
```

or

```text
lat_boundary_scheme = 2.
```

They are controlled independently at the southern and northern boundaries:

```fortran
nm1%sponge_south_width
nm1%sponge_south_timescale
nm1%sponge_north_width
nm1%sponge_north_timescale
```

Widths are in **degrees latitude inward from the actual physical model boundary**.

The sponge therefore uses the truncated domain edges, such as `65°N` and `86.5°N`, rather than incorrectly measuring distance from ±90°.

A side is disabled when either its width or its timescale is zero.

Within the sponge, `h`, `u`, and `v` are exponentially relaxed toward their initialized reference values. The damping rate uses a quadratic spatial ramp:

```text
rate = xi² / tau_boundary
```

and the update uses the exact exponential factor

```text
damp = exp(-dt * rate).
```

This makes the sponge weakest at its interior edge and strongest at the latitude boundary.

For a radiative boundary, a weak sponge can be useful for absorbing components that a simple one-dimensional characteristic boundary does not perfectly transmit, particularly slow or obliquely propagating disturbances.

Example:

```fortran
nm1%lat_boundary_scheme=2,
nm1%sponge_north_width=2.0,
nm1%sponge_north_timescale=172800.,
nm1%sponge_south_width=0.,
nm1%sponge_south_timescale=0.,
```

would apply a weak two-degree northern sponge with a two-day boundary e-folding time.

---

## 10. Nudging

The model can relax the zonal wind toward the specified reference jet:

```fortran
nm1%nudge=.true.
nm1%nudge_timescale=...
```

Two schemes are available.

### 10.1 Legacy local nudging

```fortran
nm1%nudge_scheme=0
```

nudges each longitude column independently toward `u_nudge(theta)`.

### 10.2 Global zonal-mean nudging

```fortran
nm1%nudge_scheme=1
```

nudges only the global zonal-mean zonal wind toward the reference profile.

The same zonal-mean increment is applied across all longitudes at a given latitude, so eddy departures from the zonal mean are not directly damped by the nudging term.

This mode uses a longitude-only MPI subcommunicator constructed with `MPI_CART_SUB`.

For jet-instability studies, scheme 1 is usually preferable because it maintains the background jet without directly suppressing longitudinal eddy structure.

---

## 11. Dissipation and subgrid models

Dissipation is enabled with

```fortran
nm1%viscous_dissipation=.true.
```

and selected using

```fortran
nm1%subgrid_model
```

### 11.1 Constant viscosity

```fortran
nm1%subgrid_model=1
```

uses a constant viscosity specified by

```fortran
nm1%vis
```

through the model's spherical Laplacian operator.

If

```fortran
nm1%dissipate_h=.true.
```

then layer thickness is also diffused in the constant-viscosity configuration.

### 11.2 Smagorinsky closure

```fortran
nm1%subgrid_model=2
```

uses a Smagorinsky-type dynamically varying eddy viscosity controlled by

```fortran
nm1%cvis
```

with typical values of order `0.1–0.2`.

Two Smagorinsky formulations are available.

#### Legacy scalar formulation

```fortran
nm1%smagorinsky_scheme=0
```

computes an effective viscosity and applies scalar Laplacian diffusion separately to `u` and `v`.

#### Spherical conservative SGS stress

```fortran
nm1%smagorinsky_scheme=1
```

uses a spherical, thickness-weighted SGS stress tensor and applies the divergence of that stress to the conservative momenta.

This formulation includes spherical strain/metric terms and does not directly diffuse `h`.

### 11.3 Extra equatorial meridional viscosity

The variables

```fortran
nm1%vis_eq
nm1%lat_eq
```

provide an additional viscosity contribution to `v` within a latitude band around the equator.

This is mainly a legacy/special-purpose control and is generally irrelevant for a high-latitude Saturn-only domain.

---

## 12. Planetary parameters

The main physical parameters are:

```fortran
nm1%grav
nm1%Re
nm1%rotation_period_hours
nm1%scale_height
nm1%rho
```

where:

- `grav` is gravitational acceleration in m s^-2;
- `Re` is planetary radius in metres;
- `rotation_period_hours` is the rotation period in hours;
- `scale_height` is used as the shallow-water equivalent depth;
- `rho` is retained as a model parameter but is not currently active in the shallow-water tendency equations.

Despite its historical name, `scale_height` should be interpreted as the **shallow-water equivalent depth**, not necessarily a literal atmospheric density scale height.

For Saturn, the supplied example uses an equivalent depth of `60 km`. This parameter strongly affects gravity-wave speed and therefore both the physical deformation scale and the numerical timestep constraint.

---

## 13. NetCDF output

The output file is specified using

```fortran
nm1%outputfile
```

and the output interval in seconds using

```fortran
nm1%output_interval
```

NetCDF fields include:

- `time` — simulation time;
- `phi` — longitude coordinate;
- `theta` — latitude coordinate;
- `u_nudge` — target/reference zonal-wind profile;
- `f_cor` — Coriolis parameter;
- `height` — total free-surface height `h + hs`;
- `h` — active layer thickness;
- `u` — zonal velocity;
- `v` — meridional velocity;
- `vort` — diagnosed relative vorticity.

Output variables are stored as NetCDF single-precision `REAL` values even when the model is compiled internally in double precision.

The current Makefile default is

```make
VAR_TYPE = 1
```

which selects double-precision working arithmetic through the model's numerical-type module.

---

## 14. Building the model

### Dependencies

The code requires:

- a Fortran compiler;
- MPI with Fortran bindings;
- NetCDF-Fortran;
- the `osnf` source/library directory expected by the Makefile.

The default compiler wrappers are

```make
FOR  = mpif90 -c
FOR2 = mpif90
```

and NetCDF is linked with

```make
NETCDF_LIB=-lnetcdff
```

The Makefile expects the NetCDF include/library locations to be supplied via `NETCDF_FOR` and `NETCDF_C` as appropriate for the local system.

Build with:

```bash
make
```

The resulting executable is

```text
main.exe
```

To remove compiled files:

```bash
make clean
```

or recursively clean the model and `osnf` tree with:

```bash
make cleanall
```

### Note on the supplied archive

The model Makefile expects an `osnf/` directory containing supporting numerical routines such as `numerics`, interpolation, root-finding, random-number, and other utilities. If that directory is not present in a copied/archive version of the repository, the model cannot be built until the dependency is restored.

---

## 15. Running the model

The executable expects the namelist filename as its first command-line argument:

```bash
./main.exe namelist.in
```

For MPI runs, use either `mpiexec` directly:

```bash
mpiexec -n 8 ./main.exe namelist.in
```

or the supplied helper script:

```bash
./run.sh 8 namelist.in
```

If the process-count argument is omitted, `run.sh` runs one MPI process.

The script creates an output path under `/tmp/${USER}` and writes a temporary namelist.

On systems where `/tmp/${USER}` already exists, the current script may print a harmless message such as:

```text
mkdir: /tmp/<user>: File exists
```

Using

```bash
mkdir -p /tmp/${USER}
```

would suppress this message.

---

## 16. Batch experiments

`batch_runs.sh` provides a simple parameter-sweep example.

The supplied script loops over arrays of jet speeds and Smagorinsky coefficients, edits a temporary namelist with `sed`, runs the MPI model, and writes separate NetCDF files for each case.

It is intended as an example rather than a general experiment-management system. Check the `sed` search strings against the current `namelist.in` before using it, because they depend on exact text matches.

---

## 17. Example Saturn configuration

A representative current setup is:

```fortran
&run_vars
    nm1%inputfile='saturn_winds_vs_latitude_update.nc'
    nm1%outputfile='/tmp/output.nc'

    nm1%initial_winds=1
    nm1%wind_factor=1.0
    nm1%wind_shift=0.0
    nm1%wind_reduce=0.0

    nm1%add_random_height_noise=.true.
    nm1%height_noise_scheme=1
    nm1%height_noise_amplitude=100.
    nm1%height_noise_corr_length=5.e5
    nm1%initially_geostrophic=.true.

    nm1%coriolis_scheme=1
    nm1%momentum_metric_terms=1

    nm1%lat_boundary_scheme=2
    nm1%sponge_south_width=0.
    nm1%sponge_south_timescale=0.
    nm1%sponge_north_width=0.
    nm1%sponge_north_timescale=0.

    nm1%nudge=.true.
    nm1%nudge_scheme=1
    nm1%nudge_timescale=1.7e6

    nm1%viscous_dissipation=.true.
    nm1%subgrid_model=2
    nm1%smagorinsky_scheme=1
    nm1%cvis=0.15
    nm1%dissipate_h=.false.

    nm1%runtime=7776000.
    nm1%dt=30.
    nm1%output_interval=57600.

    nm1%grav=10.44
    nm1%Re=5.45e7
    nm1%rotation_period_hours=10.6564
    nm1%scale_height=60.e3

    nm1%ip=400
    nm1%jp=240
    nm1%slat=-90.
    nm1%nlat=90.
    nm1%slat_thresh=65.
    nm1%nlat_thresh=86.5
/
```

This configuration runs a full-periodic longitude domain between 65°N and 86.5°N with observed Saturn winds, correlated height perturbations, spherical metric terms, Crank–Nicolson Coriolis treatment, zonal-mean nudging, spherical Smagorinsky stress, and radiative latitude boundaries.

---

## 18. Namelist reference

| Variable | Meaning |
|---|---|
| `inputfile` | NetCDF latitude/wind input file |
| `outputfile` | NetCDF model output path |
| `ip`, `jp` | Global longitude and latitude grid sizes |
| `runtime` | Requested integration duration, s |
| `dt` | Timestep, s |
| `output_interval` | NetCDF output interval, s |
| `grav` | Gravity, m s^-2 |
| `rho` | Density parameter; currently not active in dynamics |
| `Re` | Planetary radius, m |
| `rotation_period_hours` | Planetary rotation period, h |
| `scale_height` | Shallow-water equivalent depth, m |
| `slat`, `nlat` | Requested latitude limits, degrees |
| `slat_thresh`, `nlat_thresh` | Applied latitude clipping limits, degrees |
| `initial_winds` | `1`: Saturn profile; `2`: idealised Gaussian jet |
| `u_jet` | Idealised jet peak speed, m s^-1 |
| `theta_jet` | Idealised jet centre, degrees |
| `h_jet` | Idealised jet Gaussian width, degrees |
| `wind_factor` | Multiplicative Saturn-wind scaling |
| `wind_shift` | Saturn-profile latitude shift, degrees |
| `wind_reduce` | Velocity subtracted from Saturn profile, m s^-1 |
| `add_random_height_noise` | Enable initial height perturbations |
| `height_noise_scheme` | `0`: legacy cell noise; `1`: correlated physical-scale noise |
| `height_noise_amplitude` | Correlated-noise RMS amplitude, m |
| `height_noise_corr_length` | Correlated-noise Gaussian length scale, m |
| `initially_geostrophic` | Re-diagnose balanced wind after height perturbation |
| `coriolis_scheme` | `0`: legacy; `1`: Crank–Nicolson |
| `momentum_metric_terms` | `0`: legacy/off; `1`: spherical curvature terms |
| `lat_boundary_scheme` | `0`: legacy; `1`: free-slip wall; `2`: radiative/open |
| `sponge_south_width` | Southern sponge width inward from model edge, degrees |
| `sponge_north_width` | Northern sponge width inward from model edge, degrees |
| `sponge_south_timescale` | Southern boundary sponge e-folding time, s |
| `sponge_north_timescale` | Northern boundary sponge e-folding time, s |
| `nudge` | Enable zonal-wind nudging |
| `nudge_scheme` | `0`: local; `1`: global zonal mean |
| `nudge_timescale` | Nudging timescale, s |
| `viscous_dissipation` | Enable explicit dissipation / SGS closure |
| `subgrid_model` | `1`: constant viscosity; `2`: Smagorinsky |
| `smagorinsky_scheme` | `0`: legacy scalar; `1`: spherical conservative stress |
| `vis` | Constant viscosity coefficient |
| `cvis` | Smagorinsky coefficient |
| `dissipate_h` | Diffuse `h` for constant-viscosity model |
| `vis_eq` | Extra equatorial `v` viscosity |
| `lat_eq` | Half-width of equatorial-viscosity region, degrees |
| `restart` | Declared namelist option; restart functionality is not currently implemented |

---

## 19. Known limitations and issues

### 19.1 No implemented restart path

`nm1%restart` exists in the namelist type but is not currently used by the model driver or initialisation code. Setting it does not provide a restart capability.

A future restart implementation should read `h`, `u`, `v`, simulation time, and any required reference fields from NetCDF and continue the output record/time bookkeeping consistently.

### 19.2 No automatic CFL control

The model does not currently calculate a production CFL limit or adjust `dt` automatically.

Users must choose `dt` conservatively based on both flow velocity and gravity-wave speed, especially at high latitude where zonal grid spacing decreases as `cos(theta)`.

A useful approximate diagnostic is

```text
Cx ~ dt (|u| + sqrt(g h)) / (a cosθ dλ)
Cy ~ dt (|v| + sqrt(g h)) / (a dθ).
```

### 19.3 Lax–Wendroff is not positivity preserving

The current dynamical core does not mathematically guarantee

```text
h > 0
```

for arbitrarily strong disturbances.

If the solution becomes under-resolved or violates the timestep stability limit, a negative layer depth can occur, followed by invalid divisions in the momentum update.

The preferred response is to address resolution, timestep, forcing, boundary reflection, or dissipation—not simply clip `h`, because clipping introduces an artificial mass source.

### 19.4 Radiative boundary is approximate

`lat_boundary_scheme=2` is based on linear normal shallow-water characteristics around the initialized reference state.

It is designed to reduce reflection of outgoing disturbances from an artificial truncated latitude boundary, but it is not exactly transparent to:

- strongly nonlinear waves;
- waves incident obliquely on the boundary;
- slow balanced/Rossby disturbances;
- structures whose background state has evolved far from the initial reference profile.

A weak sponge can be combined with the radiative boundary when residual reflection matters.

### 19.5 Free-slip wall is generally artificial for truncated planetary domains

`lat_boundary_scheme=1` represents a genuine reflecting wall and can trap wave energy. It is useful as a numerical experiment but should not normally represent a latitude such as 86.5°N in Saturn's atmosphere.

For that application, scheme 2 is the more physically defensible default.

### 19.6 Latitude–longitude pole singularity

The current grid is not designed to include the geographic poles directly. Approaching ±90° makes zonal grid spacing very small and introduces coordinate singularities through `cos(theta)` and `tan(theta)`.

A truly global spherical model would benefit from dedicated polar filtering/treatment or a different grid such as a cubed sphere, Yin–Yang grid, or spectral formulation.

### 19.7 Input wind-profile assumptions

The Saturn initialisation uses `find_pos` followed by interpolation in the input latitude array. The wind input should therefore be checked for sensible ordering, duplicate latitude entries, or local reversals before extending simulations into latitude ranges not already tested.

### 19.8 `rho` currently has no dynamical effect

The density namelist parameter is stored and passed through the model but is not presently used in the shallow-water tendency calculations. It should not be interpreted as an active control on the current dynamics.

### 19.9 Runtime may extend slightly beyond the requested value

Because

```text
ntim = ceiling(runtime/dt),
```

a requested runtime that is not an integer multiple of `dt` can result in a nominal integration length slightly greater than `runtime`.

### 19.10 Batch script is text-substitution based

`batch_runs.sh` edits the namelist using literal `sed` substitutions. Changes to spacing or default values in `namelist.in` can therefore prevent a substitution from matching.

For large experiment suites, generating complete namelist files from Python or another structured tool would be more robust.

---

## 20. Recommended Saturn usage

For the current high-latitude Saturn application, a sensible starting point is:

```fortran
nm1%initial_winds=1
nm1%initially_geostrophic=.true.
nm1%momentum_metric_terms=1
nm1%coriolis_scheme=1
nm1%height_noise_scheme=1
nm1%nudge_scheme=1
nm1%lat_boundary_scheme=2
```

Start with the radiative boundary and no sponge. If outgoing disturbances still show visible reflection near the truncated northern boundary, add a weak northern sponge rather than replacing the boundary with a rigid free-slip wall.

When assessing sensitivity, vary at least:

- equivalent depth (`scale_height`);
- timestep and resolution;
- perturbation amplitude/correlation scale;
- nudging timescale;
- SGS coefficient or explicit dissipation;
- sponge width/timescale if enabled;
- position of the northern domain edge.

Boundary-condition and numerical-dissipation choices should be tested for their effect on jet morphology and wave propagation, not only whether a run remains numerically stable.

---

## 21. Source layout

The principal source files are:

```text
main.f90
```

MPI setup, namelist reading, global model setup, and call to the driver.

```text
variables.f90
```

Model, MPI, I/O, and namelist data types.

```text
initialisation.f90
```

Grid construction, Saturn/idealised wind setup, balanced height/wind initialisation, random perturbations, and metric arrays.

```text
advection.f90
```

Shallow-water dynamical core, Coriolis/metric treatment, dissipation, Smagorinsky closures, and spherical SGS stress operators.

```text
driver_code.f90
```

Time loop, output, nudging, dissipation, latitude boundary conditions, sponge layers, diagnostics, and NetCDF writing.

```text
mpi_module.f90
```

MPI datatypes, packed halo exchange, and synchronization helpers.

```text
namelist.in
```

Example Saturn configuration.

```text
run.sh
batch_runs.sh
```

Simple single-run and batch-run helper scripts.

```text
saturn_winds_vs_latitude_update.nc
```

Reference Saturn latitude/wind dataset used by the example configuration.

```text
osnf/
```

External/support numerical routines expected by the Makefile.

---

## 22. Development notes

When changing the model, useful regression checks include:

1. run the initialized balanced jet with perturbations disabled and confirm that spurious meridional motion remains small;
2. compare results across different MPI process decompositions;
3. verify global mass conservation for configurations where the latitude boundaries should be closed;
4. for radiative boundaries, monitor mass and momentum flux through the open edges rather than expecting strict closed-domain conservation;
5. compare scheme 0, scheme 1, and scheme 2 only with awareness that they represent physically different boundary behaviours;
6. test timestep and resolution sensitivity before interpreting small-scale structures as physical;
7. inspect `h_min`, maximum velocity, and boundary-region fields when diagnosing a numerical failure;
8. verify that results are not strongly dependent on sponge width/timescale when using a sponge as an absorbing layer.

For the Saturn polar-jet problem, physical realism should be judged using both numerical stability and the morphology, propagation, phase speed, and persistence of the simulated jet/wave structures.
