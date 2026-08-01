% EDDY_KINEMATICS_DRIFTERS  Vorticity, divergence and strain of the Peleliu tip
%                           vortex from a cluster of four GPS wave-buoy drifters.
%
% SCIENCE CONTEXT
%   ARCTERX-2023 Wake experiment, Peleliu (Palau) tip vortex shedding.  Four
%   CORDC Miniature Wave Buoys ("mwb" drifters) were seeded into a vortex at the
%   Peleliu tip.  The vortex detached from the island and traveled away, and
%   between 2023-05-22 05:10:54 and 2023-05-23 06:34:00 UTC all four drifters
%   traveled with it.  The goal of this script is the VORTICITY of that vortex
%   and how it evolved.
%
% WHY THIS IS DIFFERENT FROM explore_drifter_paths.m
%   explore_drifter_paths.m models each drifter's POSITION as
%
%       x_earth = x_center + radius_drifter * theta_drifter
%
%   and fits a shared center plus one radius per drifter.  Vorticity then has to
%   be inferred indirectly from the fitted radius and orbital speed, and the
%   model has a sign degeneracy: the center may sit on either side of the
%   velocity vector, so a fitted radius can come out negative.
%
%   Vorticity is a property of the VELOCITY field, and these buoys measure
%   velocity directly (GPS Doppler, variables u and v).  So here we estimate it
%   from the velocities, two independent ways, and cross-check them.
%
%   Incidentally, that older fit is LINEAR in (center, radii) once theta is
%   taken from the velocity, so fminsearch was never needed -- backslash solves
%   it in closed form.
%
% THE PHYSICS, IN PLAIN TERMS  (read this first if you are new to it)
%   Imagine dropping a tiny paddle-wheel into the water.  VORTICITY is twice the
%   rate that paddle-wheel spins.  Units are 1/s.  For water rotating like a
%   solid disc at angular rate omega, the vorticity is zeta = 2*omega.
%
%   Two ways to measure it from four drifters:
%
%   (1) CIRCULATION / STOKES' THEOREM.  Walk around the closed loop joining the
%       four drifters and add up the component of velocity pointing along your
%       direction of travel.  That total is the circulation Gamma.  Stokes'
%       theorem says Gamma equals vorticity integrated over the enclosed area:
%
%           mean vorticity inside the loop  =  Gamma / Area
%
%       Exact for the AREA-AVERAGED vorticity no matter how complicated the flow
%       inside the loop is.  That robustness is why it is the primary estimator.
%       The same trick with the OUTWARD component gives the divergence.
%
%   (2) LEAST-SQUARES VELOCITY GRADIENT.  Assume the velocity varies linearly
%       across a patch as small as the cluster:
%
%           u(x,y) = U0 + dudx*(x-xbar) + dudy*(y-ybar)
%           v(x,y) = V0 + dvdx*(x-xbar) + dvdy*(y-ybar)
%
%       Each drifter at each time is one equation, so the gradients come from a
%       linear least-squares solve -- no starting guess, no iteration, no local
%       minima.  This is the classical drifter-cluster method of Okubo and
%       Ebbesmeyer (1976), whose Eqs. (1)-(2) are exactly the model fitted in
%       fitGradient below, and independently Molinari and Kirwan (1975).
%       From the gradients:
%
%           vorticity      zeta  = dvdx - dudy      (spin)
%           divergence     delta = dudx + dvdy      (spreading out)
%           normal strain  sig_n = dudx - dvdy      (stretch along an axis)
%           shear strain   sig_s = dvdx + dudy      (shear)
%
%       and the Okubo-Weiss parameter (Okubo 1970; Weiss 1991)
%
%           OW = sig_n^2 + sig_s^2 - zeta^2
%
%       OW < 0 means rotation beats strain: you are inside a coherent vortex.
%
%   Method (1) assumes nothing about the flow inside the loop; method (2)
%   assumes it is linear across the cluster.  For a strictly linear flow the two
%   agree exactly, so ANY DISAGREEMENT MEASURES HOW NON-LINEAR THE FLOW IS
%   ACROSS THE CLUSTER.  That is information, not a nuisance -- do not average
%   the two together to "split the difference".
%
%   Molinari and Kirwan saw this in 1975: their time series were "ragged with
%   frequent changes in sign" wherever shear rates were small compared with
%   observational error.  Same failure mode, same fix -- quote the resampling
%   error bar, not the formal one.
%
% HOW WE KNOW THE ANSWER IS NOT NOISE  (the checks that matter)
%   The instinct is to compute the number and report it.  Resist that.  Every
%   number below is checked against something built from different information:
%
%   * THREE vorticity estimators: circulation, least squares, and 2 x the
%     rotation rate of the constellation (POSITIONS ONLY, no velocities).
%   * THREE divergence estimators: contour flux, least squares, and d(lnA)/dt
%     from the cluster area alone (again positions only).  For a coherent eddy
%     divergence should be near zero; it is, at ~2% of |zeta|.
%   * A LEAVE-ONE-DRIFTER-OUT (jackknife) uncertainty.  Four drifters give four
%     triangles; for a linear flow all four would agree exactly, so their spread
%     measures noise PLUS flow curvature without assuming either (Efron and
%     Gong 1983).  It comes out ~17% of |zeta|, about three times the formal
%     least-squares error -- the formal error cannot know the flow is curved.
%   * BOOTSTRAP confidence intervals on the vortex-structure fit, because a
%     single global fit is far more sensitive to binning than the medians are
%     (see the CAVEATS).
%   * AN INVARIANCE TEST on the pre-filtering.
%
% PREPROCESSING: WHY WE AVERAGE BEFORE ANYTHING ELSE
%   The buoys sample at 2 Hz and every variable -- lat, lon, u, v, height -- is
%   contaminated by surface waves.  Wave energy sits below ~30 s period and
%   peaks near 4.4 s; within a burst the wave orbital velocity has a standard
%   deviation of ~0.32 m/s, comparable to the ~0.4 m/s current we are after.
%   So we MUST low-pass before estimating anything.
%
%   Data are stored as bursts of 2048 samples (the chunk the buoy uses for its
%   onboard directional wave spectrum).  A burst is internally gap-free, but
%   successive bursts are separated by 35-202 s gaps while the buoy processes
%   and telemeters.  THEREFORE: average WITHIN a burst, never across one.  A
%   smoothing window spanning a gap is built from one-sided, sparse data and can
%   be wrong by up to 0.6 m/s at the burst edges.
%
%   ncread returns each variable as (2048 x nBurst), i.e. ONE COLUMN PER BURST,
%   so reshaping a column into blocksPerBurst blocks and averaging down them
%   never mixes data across a gap.  The default of 4 gives ~256 s averages.
%   Verified on these data: the final vorticity changes by only 0.7% when
%   blocksPerBurst runs from 1 to 16 (1024 s down to 64 s averaging).  The
%   jackknife spread does NOT shrink with longer averaging either, which proves
%   the residual uncertainty is real spatial structure, not wave noise.
%
% MODERN MATLAB USED HERE  (R2026a, full toolbox suite assumed)
%   This is a deliberate rewrite of a hand-rolled first version.  Each
%   replacement removes code that could go wrong rather than merely shortening
%   it:
%
%     timetable + retime          one common time base for four unsynchronized
%                                 drifters, replacing a manual interp1 loop
%     dictionary                  per-drifter lookup by name
%     geodetic2enu +              (Mapping Toolbox) exact WGS-84 local tangent
%       wgs84Ellipsoid            plane, replacing hand-coded radii of curvature
%     polyshape                   area and perimeter of the drifter polygon,
%                                 with self-intersection detected for us
%     fitlm                       (Statistics Toolbox) the velocity-gradient fit
%                                 WITH a coefficient covariance matrix, so the
%                                 error propagation is exact rather than assumed
%     jackknife                   (Statistics Toolbox) the built-in leave-one-out
%     fit / fittype               (Curve Fitting Toolbox) Lamb-Oseen by proper
%                                 nonlinear least squares, not a grid search
%     bootstrp + prctile          confidence intervals on that fit
%     findgroups + splitapply     radial binning without an index loop
%     smoothdata(SamplePoints=)   gap-aware moving averages on real timestamps
%     geoplot / geobasemap        the track on an actual map
%     xregion, fontsize,          modern graphics
%       colororder, exportgraphics
%
%   Speed note: fitlm is chosen over a raw backslash because we USE the
%   covariance it returns.  If you ever push this to thousands of windows, swap
%   in pagemldivide for a batched solve -- roughly 100x faster, but it gives you
%   no statistics, so you would have to propagate the errors yourself.
%
% INPUT
%   <dataDir>/mwb*_gps_timeseries.nc -- CF-1.13 trajectory files written by
%   mwb_nc.py from the raw instrument .dat bursts.
%     time    milliseconds since 1970-01-01 (UTC)
%     u,v,w   m/s, eastward/northward/upward GPS velocity
%     lat,lon degrees north / east (WGS-84)
%
% OUTPUT  (written to cfg.outDir)
%   eddy_kinematics.mat     results as timetables and a structure
%   vorticity.png           three vorticity estimators + cluster area
%   kinematics.png          divergence, strain, Okubo-Weiss
%   tracks_center.png       Earth-frame tracks, eddy frame, and a map
%   radial_profile.png      azimuthal velocity and zeta versus radius
%   eddy_frame_scatter.png  each drifter in the eddy frame, colored by time
%
% REQUIREMENTS
%   MATLAB R2026a.  Mapping, Statistics and Machine Learning, and Curve Fitting
%   Toolboxes.
%
% RESULTS FROM THESE DATA  (so you can tell whether your run is sane)
%   median zeta       -1.19e-3 1/s        (anticyclonic, i.e. clockwise)
%   Rossby number     zeta/f = -67        (f = 1.78e-5 1/s at 6.99 N)
%   rotation period   2.9 h
%   |divergence|      < 2% of |zeta|
%   Okubo-Weiss < 0   in 100% of windows -- coherent for the whole 25.4 h
%   Lamb-Oseen core   R ~ 1.1-1.2 km (quote the CI, not four digits)
%
% CAVEATS  (read before believing any of this)
%   * The cluster (~640 m across) is comparable to the vortex core (~1.2 km), so
%     the "linear flow across the cluster" assumption of method (2) is marginal.
%     That is exactly why method (1) is primary.
%   * Each vorticity value is an AREA AVERAGE over whatever the polygon spanned
%     at that moment.  The vortex is NOT solid-body -- zeta falls off with
%     radius -- so part of the apparent time variation is the cluster sampling
%     different radii.  The radial profile separates the two.
%   * THE LAMB-OSEEN FIT IS THE SENSITIVE PART.  A hand-rolled first port of
%     this script returned R = 1040 m where the reference Python implementation
%     returned 1206 m, from nothing more than different conventions for
%     percentile, moving-average edges, and window placement.  Repeating the
%     Python fit over a grid of bin counts and radius cutoffs spans
%     R = 1161-1286 m.  The medians are robust; this single global fit is not.
%     That is why R and Gamma are now reported WITH bootstrap confidence
%     intervals, and why you should quote the interval, never the point
%     estimate to four digits.
%   * The profile has not clearly turned over by the largest sampled radius
%     (~1900 m), so R is better constrained than Gamma, which extrapolates.
%   * The eddy CENTER cannot be found from an instantaneous velocity field
%     alone: a uniform background flow and a shift of the center are exactly
%     degenerate (see translationVelocity).  We break the tie by low-passing the
%     drifter-mean velocity over several orbital periods, so every center
%     estimate is conditional on that.
%   * A near-collinear cluster makes Gamma/Area blow up; epochs whose polygon
%     quality 4*pi*A/P^2 falls below cfg.minQuality are rejected.
%
% REFERENCES  (PDFs and BibTeX in ./papers, see papers/README.md)
%   Okubo, A., and C. C. Ebbesmeyer, 1976: Determination of vorticity,
%     divergence, and deformation rates from analysis of drogue observations.
%     Deep-Sea Research, 23, 349-352.  doi:10.1016/0011-7471(76)90875-5
%   Molinari, R., and A. D. Kirwan, Jr., 1975: Calculations of differential
%     kinematic properties from Lagrangian observations in the western Caribbean
%     Sea.  J. Phys. Oceanogr., 5, 483-491.
%     doi:10.1175/1520-0485(1975)005<0483:CODKPF>2.0.CO;2
%   Okubo, A., 1970: Horizontal dispersion of floatable particles in the
%     vicinity of velocity singularities such as convergences.  Deep-Sea
%     Research, 17, 445-454.  doi:10.1016/0011-7471(70)90059-8
%   Weiss, J., 1991: The dynamics of enstrophy transfer in two-dimensional
%     hydrodynamics.  Physica D, 48, 273-294.  doi:10.1016/0167-2789(91)90088-Q
%   Efron, B., and G. Gong, 1983: A leisurely look at the bootstrap, the
%     jackknife, and cross-validation.  Am. Statistician, 37(1), 36-48.
%     doi:10.1080/00031305.1983.10483087
%   Saffman, P. G., 1992: Vortex Dynamics.  Cambridge University Press.

%% ------------------------------------------------------------------------
%  Configuration
%  ------------------------------------------------------------------------

cfg = struct(); % Configuration parameters
cfg.dataDir        = string(fileparts(mfilename("fullpath")));
cfg.outDir         = fullfile(cfg.dataDir, "eddy_out_matlab");
cfg.drifters       = append("mwb", ["458d02", "788d01", "790d01", "793d02"]);
cfg.tStart         = datetime(2023,5,22,05,10,54, TimeZone="UTC");
cfg.tEnd           = datetime(2023,5,23,06,34,00, TimeZone="UTC");
cfg.blocksPerBurst = 4;             % 2048/4 = 512 samples = 256 s averages
cfg.window         = minutes(30);   % least-squares window width
cfg.step           = minutes(10);   % spacing between successive windows
cfg.transPeriod    = hours(6);      % low-pass for eddy translation (~2 orbits)
cfg.minQuality     = 0.10;          % reject polygons flatter than 4*pi*A/P^2
cfg.maxGap         = minutes(30);   % never interpolate a drifter this far
cfg.snrCenter      = 3;             % require |zeta| > snr * its formal error
cfg.maxDispScales  = 3;             % reject a center this many cluster widths out
cfg.nRadialBins    = 14;
cfg.minPerBin      = 20;            % ignore radial bins sparser than this
cfg.nBoot          = 1000;          % bootstrap replicates for the Oseen CI

if ~isfolder(cfg.outDir); mkdir(cfg.outDir); end
rng(20230522, "twister");           % reproducible bootstrap

% polyshape repairs duplicate and collinear vertices routinely here; that is
% wanted behavior, not a problem, so silence the per-call chatter.
warnState = warning("off", "MATLAB:polyshape:repairedBySimplify");
cleanupWarn = onCleanup(@() warning(warnState));

%% ------------------------------------------------------------------------
%  Load, block-average, and put every drifter on one common time base
%  ------------------------------------------------------------------------
%  The buoys are NOT synchronized -- each runs its own sample/process/telemeter
%  cycle, so their bursts drift out of phase.  Vorticity needs all four
%  positions and velocities at the SAME instant.  Holding each drifter in a
%  timetable makes retime the natural one-line answer to that.

nD = numel(cfg.drifters);
perDrifter = dictionary();
for name = cfg.drifters
    fn = fullfile(cfg.dataDir, name + "_gps_timeseries.nc");
    perDrifter(name) = {loadBlocks(fn, cfg.tStart, cfg.tEnd, cfg.blocksPerBurst)};
end % for name

% Common regular grid spanning the interval all four have in common
starts = NaT(1,nD, TimeZone="UTC");  stops = starts;  steps = seconds(nan(1,nD));
for k = 1:nD
    tt = perDrifter{cfg.drifters(k)};
    starts(k) = tt.Time(1);  stops(k) = tt.Time(end);
    steps(k)  = median(diff(tt.Time));
end % for k
step  = median(steps);
tGrid = (max(starts) : step : min(stops)).';
nT    = numel(tGrid);

[lat, lon, U, V] = deal(nan(nT, nD));
for k = 1:nD
    tt = perDrifter{cfg.drifters(k)};
    rt = retime(tt, tGrid, "linear");
    % Mark, rather than silently accept, any grid point that had to reach too
    % far for a real observation of this drifter.
    far = min(abs(tGrid - tt.Time.'), [], 2) > cfg.maxGap;
    rt{far, :} = NaN;
    lat(:,k) = rt.lat;  lon(:,k) = rt.lon;
    U(:,k) = rt.u;      V(:,k) = rt.v;
end % for k

% One local tangent plane for the whole cluster.  All drifters MUST share the
% origin, because vorticity depends only on the relative geometry.
lat0 = mean(lat, "all", "omitnan");
lon0 = mean(lon, "all", "omitnan");
[X, Y] = geodetic2enu(lat, lon, 0, lat0, lon0, 0, wgs84Ellipsoid("meter"));
f = coriolis(lat0);

% RowTimes= (rather than passing the times positionally) guarantees the time
% dimension is named "Time".  Passing them first names it after the input
% variable instead -- "tGrid" -- and every later TT.Time then errors.
TT = timetable(X, Y, U, V, lat, lon, RowTimes=tGrid, ...
    VariableNames=["X" "Y" "U" "V" "lat" "lon"]);

fprintf("Common grid: %d epochs at %s, %s to %s\n", ...
    nT, string(step), string(tGrid(1)), string(tGrid(end)));

%% ------------------------------------------------------------------------
%  Estimator 1 (PRIMARY): circulation and flux around the drifter polygon
%  ------------------------------------------------------------------------

circ = circulationKinematics(X, Y, U, V, cfg.minQuality);
jack = circulationJackknife(X, Y, U, V, cfg.minQuality);

%% ------------------------------------------------------------------------
%  Estimator 2: least-squares velocity gradient
%  ------------------------------------------------------------------------

fitTT = fitGradient(tGrid, X, Y, U, V, cfg.window, cfg.step);

%% ------------------------------------------------------------------------
%  Estimator 3: constellation rotation, from POSITIONS ONLY
%  ------------------------------------------------------------------------

rot = rotationRate(tGrid, X, Y);

%% ------------------------------------------------------------------------
%  Independent divergence check from the cluster area alone
%  ------------------------------------------------------------------------

divArea = divergenceFromArea(tGrid, circ.area, hours(1));

%% ------------------------------------------------------------------------
%  Eddy translation and center
%  ------------------------------------------------------------------------

trans  = translationVelocity(tGrid, U, V, cfg.transPeriod);
scaleC = interp1(tGrid, sqrt(abs(circ.area)), fitTT.Time);
[cx, cy] = eddyCenter(fitTT, trans, tGrid, cfg.snrCenter, ...
    cfg.maxDispScales, scaleC);
fitTT.center_x = cx;
fitTT.center_y = cy;

fprintf("Eddy center: %d/%d windows pass the quality gates (%.0f%%)\n", ...
    sum(isfinite(cx)), numel(cx), 100*mean(isfinite(cx)));

%% ------------------------------------------------------------------------
%  Radial structure of the vortex
%  ------------------------------------------------------------------------
%  Each circulation value above mixes radii, being a mean over the whole
%  polygon.  Referred to the fitted center, every drifter sample instead becomes
%  one (radius, azimuthal velocity) pair, and pooling them resolves the
%  structure the area-average hides.

good = isfinite(cx);
cxi  = interp1(fitTT.Time(good), cx(good), tGrid, "linear", NaN);
cyi  = interp1(fitTT.Time(good), cy(good), tGrid, "linear", NaN);
prof = radialProfile(X, Y, U, V, cxi, cyi, trans, cfg.nRadialBins, cfg.minPerBin);
oseen = fitOseen(prof, cfg.nBoot, cfg.minPerBin);

%% ------------------------------------------------------------------------
%  Report
%  ------------------------------------------------------------------------

zFinite = circ.zeta(isfinite(circ.zeta));
zMed = median(zFinite);

fprintf("\n--- VORTICITY (circulation / Stokes, primary) ---\n");
fprintf("  median            : %+.3e 1/s   Rossby %+.1f\n", zMed, zMed/f);
fprintf("  MAD               : %.3e 1/s\n", median(abs(zFinite - zMed)));
fprintf("  rotation period   : %.2f h  (as solid body, T = 4*pi/|zeta|)\n", ...
    hours(seconds(2*pi/abs(zMed/2))));
fprintf("  leave-one-out 1sd : %.3e 1/s  (%.0f%% of |zeta|)\n", ...
    median(jack.spread,"omitnan"), ...
    100*median(jack.spread./abs(jack.zeta),"omitnan"));

fprintf("\n--- CONSISTENCY ---\n");
fprintf("  constellation turned %+.2f revolutions\n", rot.turns);
fprintf("  |divergence|/|zeta| : contour %.3f, area %.3f, lsq %.3f\n", ...
    abs(median(circ.delta,"omitnan"))/abs(zMed), ...
    abs(median(divArea,"omitnan"))/abs(zMed), ...
    abs(median(fitTT.delta,"omitnan"))/abs(zMed));
fprintf("  Okubo-Weiss < 0 in %.0f%% of windows\n", ...
    100*mean(fitTT.OW(isfinite(fitTT.OW)) < 0));

fprintf("\n--- VORTEX STRUCTURE (Lamb-Oseen, %d bootstrap replicates) ---\n", ...
    oseen.nBootOK);
fprintf("  core radius R     : %.0f m   95%% CI [%.0f, %.0f]\n", ...
    oseen.radius, oseen.radiusCI(1), oseen.radiusCI(2));
fprintf("  circulation Gamma : %.0f m^2/s   95%% CI [%.0f, %.0f]\n", ...
    oseen.gamma, oseen.gammaCI(1), oseen.gammaCI(2));
fprintf("  core vorticity    : %+.3e 1/s  (Rossby %+.0f)\n", ...
    oseen.zetaCore, oseen.zetaCore/f);
fprintf("  adjusted R-squared: %.3f\n", oseen.rsquare);

results.cfg = cfg;   results.lat0 = lat0;  results.lon0 = lon0;  results.f = f;
results.circ = circ; results.jack = jack;  results.rot = rot;
results.divArea = divArea;  results.trans = trans;
results.prof = prof; results.oseen = oseen;
results.centerX = cxi;      results.centerY = cyi;

save(fullfile(cfg.outDir, "eddy_kinematics.mat"), ...
    "TT", "fitTT", "results", "-v7.3");

%% ------------------------------------------------------------------------
%  Figures
%  ------------------------------------------------------------------------

makeFigures(tGrid, fitTT, X, Y, lat, lon, cfg, circ, jack, rot, divArea, ...
    cxi, cyi, prof, oseen, f);

fprintf("\nWrote results and figures to %s\n", cfg.outDir);

%% ========================================================================
%  Local functions
%  ========================================================================

function tt = loadBlocks(fn, tStart, tEnd, blocksPerBurst)
% LOADBLOCKS  Read one drifter NetCDF, block-average each burst, return a timetable.
%
%   ncread returns (sampling_dwell x sampling_period) = (2048 x nBurst), so each
%   COLUMN is one contiguous, gap-free burst.  Reshaping a column into
%   blocksPerBurst blocks and averaging down them therefore never mixes data
%   across an inter-burst gap -- which is the whole point (see the header).
arguments (Input)
    fn (1,1) string {mustBeFile}
    tStart (1,1) datetime
    tEnd (1,1) datetime
    blocksPerBurst (1,1) double {mustBePositive, mustBeInteger}
end % arguments Input
arguments (Output)
    tt timetable
end % arguments Output

fprintf("Reading %s\n", fn);
names = ["time", "lat", "lon", "u", "v"];
raw = struct();
for name = names
    raw.(name) = double(ncread(fn, name));      % (2048 x nBurst)
end % for name

nSamp = size(raw.time, 1);
if mod(nSamp, blocksPerBurst) ~= 0
    error("eddy:blockSize", "%d samples per burst is not divisible by %d", ...
        nSamp, blocksPerBurst);
end % if mod

% Keep whole bursts that overlap the interval.  Trimming by BURST rather than by
% sample preserves the gap-free block structure we are about to exploit.
t0 = posixtime(tStart)*1000;
t1 = posixtime(tEnd)*1000;
keep = raw.time(end,:) >= t0 & raw.time(1,:) <= t1;
if ~any(keep)
    error("eddy:noBursts", "%s: no bursts inside the requested interval", fn);
end % if ~any

m = nSamp / blocksPerBurst;
blocked = struct();
for name = names
    a = reshape(raw.(name)(:,keep), m, blocksPerBurst, []);
    blocked.(name) = reshape(mean(a, 1), [], 1);
end % for name

inside = blocked.time >= t0 & blocked.time <= t1;
% RowTimes= for the same reason as elsewhere: it pins the time dimension name
% to "Time" rather than leaving it to be inferred from the input expression.
rowTimes = datetime(blocked.time(inside)/1000, ConvertFrom="posixtime", ...
    TimeZone="UTC");
tt = timetable(blocked.lat(inside), blocked.lon(inside), ...
    blocked.u(inside), blocked.v(inside), ...
    RowTimes=rowTimes, VariableNames=["lat","lon","u","v"]);
end % loadBlocks

function f = coriolis(lat0)
% CORIOLIS  Planetary vorticity f = 2*Omega*sin(latitude), in 1/s.
%   Vorticity is only meaningful relative to something, and f is the natural
%   yardstick.  zeta/f is the Rossby number: >> 1 means the eddy's own rotation
%   utterly dominates the Earth's, which is the submesoscale regime this vortex
%   lives in.
arguments (Input)
    lat0 (1,1) double
end % arguments Input
arguments (Output)
    f (1,1) double
end % arguments Output
f = 2 * 7.292115e-5 * sind(lat0);
end % coriolis

function [q, A, P] = polygonQuality(x, y)
% POLYGONQUALITY  Isoperimetric quotient 4*pi*A/P^2, plus area and perimeter.
%   1 for a circle, ~0.785 for a square, tending to 0 as the drifters approach a
%   straight line.  Because vorticity is circulation DIVIDED BY AREA, a sliver
%   cluster divides a small noisy circulation by a vanishing area and the
%   estimate explodes.  This quantity is scale-free, so it rejects degenerate
%   SHAPES without rejecting a merely small cluster.
%
%   polyshape does the geometry -- it detects and repairs self-intersection and
%   returns a consistent area.  Doing that by hand is where sign errors live.
arguments (Input)
    x (:,1) double
    y (:,1) double
end % arguments Input
arguments (Output)
    q (1,1) double {mustBeNonnegative}    % isoperimetric quotient, 0..1
    A (1,1) double {mustBeNonnegative}    % area, m^2
    P (1,1) double {mustBeNonnegative}    % perimeter, m
end % arguments Output

[q, A, P] = deal(0);
if numel(x) < 3 || any(~isfinite([x; y])); return; end
pgon = polyshape(x, y, Simplify=true);
if pgon.NumRegions ~= 1; return; end
A = area(pgon);
P = perimeter(pgon);
if P > 0; q = 4*pi*A/P^2; end
end % polygonQuality

function out = circulationKinematics(X, Y, U, V, minQuality)
% CIRCULATIONKINEMATICS  Vorticity and divergence by contour integral.
%
%   Stokes' theorem:      zeta_mean  = (1/A) * closed integral of u . dl
%   Divergence theorem:   delta_mean = (1/A) * closed integral of u . n dl
%
%   Both use the trapezoidal rule along each polygon edge.  The only
%   approximation is that velocity varies linearly ALONG an edge; nothing is
%   assumed about the interior, which is what makes this robust when the cluster
%   is not small compared with the vortex.
arguments (Input)
    X double
    Y double
    U double
    V double
    minQuality (1,1) double
end % arguments Input
arguments (Output)
    out (1,1) struct                      % fields zeta, delta, area, quality
end % arguments Output

nT = size(X,1);
[zeta, delta, areaOut, quality] = deal(nan(nT,1));

for k = 1:nT
    ok = isfinite(X(k,:)) & isfinite(Y(k,:)) & isfinite(U(k,:)) & isfinite(V(k,:));
    if nnz(ok) < 3; continue; end
    [zeta(k), delta(k), areaOut(k), quality(k)] = ...
        loopIntegral(X(k,ok).', Y(k,ok).', U(k,ok).', V(k,ok).', minQuality);
end % for k

out = struct("zeta", zeta, "delta", delta, "area", areaOut, "quality", quality);
end % circulationKinematics

function [zeta, delta, A, q] = loopIntegral(xs, ys, us, vs, minQuality)
% LOOPINTEGRAL  One epoch of the contour integrals.  Shared with the jackknife.
arguments (Input)
    xs (:,1) double
    ys (:,1) double
    us (:,1) double
    vs (:,1) double
    minQuality (1,1) double
end % arguments Input
arguments (Output)
    zeta (1,1) double                     % vorticity, 1/s (NaN if rejected)
    delta (1,1) double                    % divergence, 1/s (NaN if rejected)
    A (1,1) double {mustBeNonnegative}    % polygon area, m^2
    q (1,1) double {mustBeNonnegative}    % polygon quality, 0..1
end % arguments Output
[zeta, delta] = deal(NaN);
[A, q] = deal(0);

% Order the vertices counterclockwise about their centroid.  Out of order, the
% polygon self-intersects and both area and circulation come out wrong.
[~, order] = sort(atan2(ys - mean(ys), xs - mean(xs)));
xs = xs(order); ys = ys(order); us = us(order); vs = vs(order);

[q, A] = polygonQuality(xs, ys);
if A < 1 || q < minQuality; return; end

xn = circshift(xs,-1); yn = circshift(ys,-1);
un = circshift(us,-1); vn = circshift(vs,-1);
dx = xn - xs;  dy = yn - ys;

circulation = sum(0.5*(us+un).*dx + 0.5*(vs+vn).*dy);   % tangential
flux        = sum(0.5*(us+un).*dy - 0.5*(vs+vn).*dx);   % outward normal
zeta  = circulation / A;
delta = flux / A;
end % loopIntegral

function out = circulationJackknife(X, Y, U, V, minQuality)
% CIRCULATIONJACKKNIFE  Leave-one-drifter-out spread of the vorticity.
%
%   Four drifters give four distinct triangles.  If the velocity field really
%   were linear across the cluster, every triangle would return exactly the same
%   vorticity, so their spread measures measurement noise PLUS curvature of the
%   flow -- without assuming a noise model or linearity (Efron and Gong 1983).
%   THIS is the error bar to quote: it is about three times the formal
%   least-squares error, which cannot know the flow is curved.
%
%   The Statistics Toolbox jackknife() resamples by ROWS, so each epoch is
%   packed as an (nDrifters x 4) matrix of [x y u v].  Using the built-in rather
%   than an index loop is not merely shorter: it makes the leave-one-out
%   structure explicit to the reader.
arguments (Input)
    X double
    Y double
    U double
    V double
    minQuality (1,1) double
end % arguments Input
arguments (Output)
    out (1,1) struct                      % fields all, zeta, spread
end % arguments Output

[nT, nD] = size(X);
all4 = nan(nT, nD);

for k = 1:nT
    M = [X(k,:).', Y(k,:).', U(k,:).', V(k,:).'];
    if any(~isfinite(M), "all"); continue; end
    all4(k,:) = jackknife(@(m) ...
        loopIntegral(m(:,1), m(:,2), m(:,3), m(:,4), minQuality), M).';
end % for k

% Median and MAD, not mean and standard deviation: with only four triangles a
% single near-degenerate one would otherwise dominate the answer.
med = median(all4, 2, "omitnan");
out = struct("all", all4, "zeta", med, ...
    "spread", 1.4826*median(abs(all4 - med), 2, "omitnan"));
end % circulationJackknife

function fitTT = fitGradient(t, X, Y, U, V, window, step)
% FITGRADIENT  Sliding-window least-squares fit of the affine velocity field.
%
%   For each window we fit, by least squares,
%
%       u = U0 + dudx*dx + dudy*dy + dUdt*dt
%       v = V0 + dvdx*dx + dvdy*dy + dVdt*dt
%
%   The dt column lets the background flow drift linearly through the window
%   instead of forcing that signal into the gradients.
%
%   fitlm rather than backslash, because we USE what it returns: the coefficient
%   covariance gives an exact variance for zeta = dvdx - dudy rather than an
%   assumed one.  (For very many windows, pagemldivide does a batched solve
%   roughly 100x faster but returns no statistics.)
%
%   Windows are placed so every one is FULL -- no partial windows at the ends.
%   That is why this returns slightly fewer windows than an implementation that
%   lets the first and last run short on data.
arguments (Input)
    t datetime
    X double
    Y double
    U double
    V double
    window (1,1) duration
    step (1,1) duration
end % arguments Input
arguments (Output)
    fitTT timetable
end % arguments Output

centers = ((t(1) + window/2) : step : (t(end) - window/2)).';
nW = numel(centers);
[zeta, delta, sigN, sigS, OW, U0, V0, zetaErr, deltaErr, rmse, xc, yc, nObs] ...
    = deal(nan(nW,1));
grad = nan(nW,2,2);

for i = 1:nW
    inWin = abs(t - centers(i)) <= window/2;
    dt = seconds(t(inWin) - centers(i));
    xx = X(inWin,:); yy = Y(inWin,:); uu = U(inWin,:); vv = V(inWin,:);
    tt = repmat(dt, 1, size(xx,2));

    ok = isfinite(xx) & isfinite(yy) & isfinite(uu) & isfinite(vv);
    if nnz(any(ok,1)) < 3 || nnz(ok) < 6; continue; end

    xs = xx(ok); ys = yy(ok);
    xbar = mean(xs); ybar = mean(ys);
    D = [xs - xbar, ys - ybar, tt(ok)];

    mdlU = fitlm(D, uu(ok), VarNames=["dx","dy","dt","u"]);
    mdlV = fitlm(D, vv(ok), VarNames=["dx","dy","dt","v"]);

    cU = mdlU.Coefficients.Estimate;    % [intercept; dudx; dudy; dUdt]
    cV = mdlV.Coefficients.Estimate;
    CU = mdlU.CoefficientCovariance;
    CV = mdlV.CoefficientCovariance;

    dudx = cU(2); dudy = cU(3);
    dvdx = cV(2); dvdy = cV(3);

    zeta(i)  = dvdx - dudy;
    delta(i) = dudx + dvdy;
    sigN(i)  = dudx - dvdy;
    sigS(i)  = dvdx + dudy;
    OW(i)    = sigN(i)^2 + sigS(i)^2 - zeta(i)^2;

    % u and v are fitted independently, so their variances add.
    zetaErr(i)  = sqrt(CV(2,2) + CU(3,3));
    deltaErr(i) = sqrt(CU(2,2) + CV(3,3));

    U0(i) = cU(1); V0(i) = cV(1);
    rmse(i) = hypot(mdlU.RMSE, mdlV.RMSE)/sqrt(2);
    xc(i) = xbar; yc(i) = ybar; nObs(i) = nnz(ok);
    grad(i,:,:) = [dudx dudy; dvdx dvdy];
end % for i

% RowTimes= keeps the time dimension named "Time"; passing centers first would
% name it "centers" and break every fitTT.Time reference downstream.
fitTT = timetable(zeta, delta, sigN, sigS, OW, U0, V0, ...
    zetaErr, deltaErr, rmse, xc, yc, nObs, RowTimes=centers);
fitTT = addprop(fitTT, "grad", "table");
fitTT.Properties.CustomProperties.grad = grad;
end % fitGradient

function out = rotationRate(t, X, Y)
% ROTATIONRATE  Turning of the drifter constellation, from POSITIONS ONLY.
%
%   Tracks the mean angle of the drifters about their own centroid, unwrapped so
%   full turns accumulate.  For solid-body rotation d(theta)/dt = zeta/2, so
%   2*omega is directly comparable with the other two estimators -- and it uses
%   no velocity data at all, which makes it genuinely independent.
%
%   On these data the constellation turned -6.86 revolutions where zeta/2
%   integrated predicts -8.95.  That 23% shortfall is NOT an error: it is the
%   signature of the vortex not being solid-body.  Drifters further out orbit
%   more slowly than the core vorticity implies.
arguments (Input)
    t datetime
    X double
    Y double
end % arguments Input
arguments (Output)
    out (1,1) struct                      % fields theta, omega, turns
end % arguments Output

xc = mean(X, 2, "omitnan");
yc = mean(Y, 2, "omitnan");
theta = mean(unwrap(atan2(Y - yc, X - xc), [], 1), 2, "omitnan");
ts = seconds(t - t(1));
out = struct("theta", theta, "omega", gradient(theta, ts), ...
    "turns", (theta(end) - theta(1))/(2*pi));
end % rotationRate

function out = divergenceFromArea(t, areaIn, smoothWidth)
% DIVERGENCEFROMAREA  delta = d(ln A)/dt, from POSITIONS ONLY.
%
%   A patch of fluid spreading out has positive divergence, and its area grows
%   at exactly d(ln A)/dt.  Using only positions makes this independent of both
%   other divergence estimates.  All three agreeing that divergence is a couple
%   of percent of |zeta| is what lets us call the flow effectively
%   two-dimensional and non-divergent.
%
%   smoothdata with SamplePoints respects the actual timestamps, so an irregular
%   or gappy series is handled correctly instead of being treated as evenly
%   spaced.
arguments (Input)
    t datetime
    areaIn (:,1) double
    smoothWidth (1,1) duration
end % arguments Input
arguments (Output)
    out (:,1) double                      % divergence, 1/s
end % arguments Output

out = nan(size(areaIn));
ok = isfinite(areaIn) & areaIn > 0;
if nnz(ok) < 5; return; end
la = smoothdata(log(areaIn(ok)), "movmean", smoothWidth, SamplePoints=t(ok));
out(ok) = gradient(la, seconds(t(ok) - t(1)));
end % divergenceFromArea

function trans = translationVelocity(t, U, V, period)
% TRANSLATIONVELOCITY  Eddy translation, as the heavily low-passed cluster mean.
%
%   IMPORTANT CONCEPTUAL POINT.  At a single instant you CANNOT separate the
%   position of the vortex center from a uniform background current.  Write the
%   solid-body field about a center c:
%
%       u = omega * zhat x (x - c) = omega * zhat x x  -  omega * zhat x c
%
%   The second term is a CONSTANT vector, indistinguishable from a uniform
%   translation.  Moving the center and adding background flow are the same
%   thing to the data.  The degeneracy is exact, not merely ill-conditioned.
%
%   We break it with time: averaging over several orbital periods cancels the
%   rotational part (it goes round and sums to zero) and leaves the translation.
%   Every center estimate is CONDITIONAL on this -- if the center track looks
%   wrong, suspect this assumption first.
arguments (Input)
    t datetime
    U double
    V double
    period (1,1) duration
end % arguments Input
arguments (Output)
    trans (:,2) double                    % [east north] translation, m/s
end % arguments Output

mu = mean(U, 2, "omitnan");
mv = mean(V, 2, "omitnan");
trans = [smoothdata(mu, "movmean", period, "omitmissing", SamplePoints=t), ...
         smoothdata(mv, "movmean", period, "omitmissing", SamplePoints=t)];
end % translationVelocity

function [cx, cy] = eddyCenter(fitTT, trans, t, snr, maxScales, clusterScale)
% EDDYCENTER  The elliptic critical point of each windowed affine fit.
%
%   In a frame moving with the eddy the center is where the velocity vanishes:
%
%       A * (x - xbar) = c - U0     =>     x = xbar + A \ (c - U0)
%
%   Three gates, because A\(...) is a division and divisions misbehave:
%     1. Okubo-Weiss < 0, so the critical point is a center (spiral) and not a
%        saddle.  A saddle is a fine critical point but it is not the middle of
%        a vortex.
%     2. |zeta| resolved at snr sigma.  A^-1 scales as 1/|zeta|, so where the
%        rotation is weak a small velocity mismatch becomes an enormous
%        displacement.  Omitting this gate produced 13 km center excursions in
%        testing.
%     3. Displacement within maxScales cluster widths.  Drifters all sitting on
%        one side of a distant point do not constrain it.
arguments (Input)
    fitTT timetable
    trans double
    t datetime
    snr (1,1) double
    maxScales (1,1) double
    clusterScale (:,1) double
end % arguments Input
arguments (Output)
    cx (:,1) double                       % center east of origin, m
    cy (:,1) double                       % center north of origin, m
end % arguments Output

n = height(fitTT);
[cx, cy] = deal(nan(n,1));
grad = fitTT.Properties.CustomProperties.grad;
ct = interp1(t, trans, fitTT.Time);

for i = 1:n
    A = squeeze(grad(i,:,:));
    if ~all(isfinite(A), "all") || ~(fitTT.OW(i) < 0); continue; end
    if ~isfinite(fitTT.zetaErr(i)) || abs(fitTT.zeta(i)) < snr*fitTT.zetaErr(i)
        continue
    end % if ~isfinite
    if abs(det(A)) < 1e-14; continue; end

    d = A \ [ct(i,1) - fitTT.U0(i); ct(i,2) - fitTT.V0(i)];
    if isfinite(clusterScale(i)) && hypot(d(1), d(2)) > maxScales*clusterScale(i)
        continue
    end % if isfinite
    cx(i) = fitTT.xc(i) + d(1);
    cy(i) = fitTT.yc(i) + d(2);
end % for i
end % eddyCenter

function prof = radialProfile(X, Y, U, V, cx, cy, trans, nBin, minPerBin)
% RADIALPROFILE  Azimuthal velocity and local vorticity versus radius.
%
%   Referred to the fitted center, every drifter sample is one (r, v_theta)
%   pair.  Pooling them over the record resolves the radial structure that the
%   polygon area-average necessarily smears out:
%
%       Gamma(r) = 2*pi*r*v_theta(r)         circulation within radius r
%       zeta(r)  = (1/r) d(r v_theta)/dr     local vorticity
%
%   Radial velocity v_r is returned too.  For a coherent, non-dispersing vortex
%   it should be small compared with v_theta.  That is a CHECK, not an
%   assumption -- if v_r came out comparable to v_theta, the center would be
%   wrong or the vortex would be falling apart.
%
%   discretize + findgroups + splitapply do the binning, which keeps the edges
%   explicit and the per-bin statistics inspectable afterwards.
arguments (Input)
    X double
    Y double
    U double
    V double
    cx (:,1) double
    cy (:,1) double
    trans double
    nBin (1,1) double {mustBePositive, mustBeInteger}
    minPerBin (1,1) double {mustBePositive, mustBeInteger}
end % arguments Input
arguments (Output)
    prof struct
end % arguments Output

dx = X - cx;  dy = Y - cy;
r  = hypot(dx, dy);
th = atan2(dy, dx);
du = U - trans(:,1);  dv = V - trans(:,2);
vt = -du.*sin(th) + dv.*cos(th);      % azimuthal, positive counterclockwise
vr =  du.*cos(th) + dv.*sin(th);      % radial, positive outward

ok = isfinite(r) & isfinite(vt) & isfinite(vr);
T = table(r(ok), vt(ok), vr(ok), VariableNames=["r","vt","vr"]);

edges = linspace(0, prctile(T.r, 97), nBin+1);
T.bin = discretize(T.r, edges);
T = rmmissing(T, DataVariables="bin");

[G, binID] = findgroups(T.bin);
count = splitapply(@numel, T.vt, G);
keep  = count >= minPerBin;

prof.r         = splitapply(@median, T.r,  G);
prof.vTheta    = splitapply(@median, T.vt, G);
prof.vR        = splitapply(@median, T.vr, G);
prof.vThetaErr = splitapply(@(x) 1.4826*median(abs(x - median(x)))/sqrt(numel(x)), ...
                            T.vt, G);
prof.n         = count;
prof.binID     = binID;

fn = ["r" "vTheta" "vR" "vThetaErr" "n" "binID"];
for name = fn
    prof.(name) = prof.(name)(keep);
end % for name

prof.edges     = edges;
prof.minPerBin = minPerBin;
prof.raw       = T;                   % kept so the bootstrap can resample it

prof.zeta = nan(size(prof.r));
if numel(prof.r) > 2
    prof.zeta = gradient(prof.r .* prof.vTheta, prof.r) ./ prof.r;
end % if numel
end % radialProfile

function out = fitOseen(prof, nBoot, minPerBin)
% FITOSEEN  Lamb-Oseen vortex fitted to the azimuthal velocity profile.
%
%       v_theta(r) = Gamma/(2*pi*r) * (1 - exp(-r^2/R^2))
%
%   The classical viscously-spreading line vortex (Saffman 1992): solid body
%   near the middle, falling off like 1/r outside a core of radius R.  Its
%   vorticity is a Gaussian, Gamma/(pi R^2) * exp(-r^2/R^2).
%
%   WHY THE BOOTSTRAP IS NOT OPTIONAL HERE.  Unlike the medians, this is a
%   SINGLE GLOBAL FIT to about a dozen binned points, and it is genuinely
%   sensitive to how those bins were formed.  Two correct implementations of
%   this analysis, differing only in percentile convention and moving-average
%   edge handling, returned R = 1040 m and R = 1206 m.  A point estimate quoted
%   to four digits would be false precision.  So we resample the underlying
%   (r, v_theta) pairs, RE-BIN each replicate -- so binning sensitivity is
%   inside the interval rather than conditioned away -- refit, and report
%   percentile confidence intervals.  Quote the interval.
arguments (Input)
    prof struct
    nBoot (1,1) double {mustBePositive, mustBeInteger} = 1000
    minPerBin (1,1) double {mustBePositive, mustBeInteger} = 20
end % arguments Input
arguments (Output)
    out struct
end % arguments Output

ft = fittype(@(Gamma, R, r) Gamma./(2*pi*r) .* (1 - exp(-(r./R).^2)), ...
    independent="r", coefficients={'Gamma','R'});

ok = isfinite(prof.r) & isfinite(prof.vTheta) & prof.r > 0;
opts = fitoptions(ft);
opts.StartPoint = [-6000, 1200];
opts.Lower      = [-1e7, 50];
opts.Upper      = [ 1e7, 1e5];
opts.Weights    = 1 ./ max(prof.vThetaErr(ok), 1e-9).^2;

[fobj, gof] = fit(prof.r(ok), prof.vTheta(ok), ft, opts);

out.gamma    = fobj.Gamma;
out.radius   = fobj.R;
out.zetaCore = fobj.Gamma / (pi * fobj.R^2);
out.rsquare  = gof.adjrsquare;
out.fitobj   = fobj;

% confint gives the interval from the fit's own covariance; the bootstrap below
% is wider because it also captures binning choice, which confint cannot see.
try
    ci = confint(fobj, 0.95);
    out.gammaCIfit  = ci(:,1);
    out.radiusCIfit = ci(:,2);
catch
    out.gammaCIfit  = [NaN; NaN];
    out.radiusCIfit = [NaN; NaN];
end % try

bootOpts = opts;
bootOpts.Weights = [];                % re-binned replicates carry no error bars
boot = bootstrp(nBoot, ...
    @(idx) refitOseen(prof.raw(idx,:), prof.edges, ft, bootOpts, minPerBin), ...
    (1:height(prof.raw)).');
boot = boot(all(isfinite(boot), 2), :);

out.nBootOK  = size(boot,1);
if out.nBootOK >= 20
    out.gammaCI  = prctile(boot(:,1), [2.5 97.5]).';
    out.radiusCI = prctile(boot(:,2), [2.5 97.5]).';
else
    [out.gammaCI, out.radiusCI] = deal([NaN; NaN]);
end % if out.nBootOK
out.boot = boot;
end % fitOseen

function p = refitOseen(T, edges, ft, opts, minPerBin)
% REFITOSEEN  One bootstrap replicate: re-bin the resampled pairs and refit.
arguments (Input)
    T table
    edges (1,:) double
    ft fittype
    opts                                  % fitoptions object (class name is not a documented API)
    minPerBin (1,1) double {mustBePositive, mustBeInteger}
end % arguments Input
arguments (Output)
    p (1,2) double                        % [Gamma R]; [NaN NaN] if degenerate
end % arguments Output
p = [NaN NaN];
T.bin = discretize(T.r, edges);
T = rmmissing(T, DataVariables="bin");
if height(T) < 4*minPerBin; return; end

G = findgroups(T.bin);
count = splitapply(@numel, T.vt, G);
rb = splitapply(@median, T.r,  G);
vb = splitapply(@median, T.vt, G);
keep = count >= minPerBin & isfinite(rb) & isfinite(vb) & rb > 0;
if nnz(keep) < 4; return; end

try
    fobj = fit(rb(keep), vb(keep), ft, opts);
    p = [fobj.Gamma, fobj.R];
catch
    % a degenerate replicate simply drops out of the interval
end % try
end % refitOseen

function makeFigures(t, fitTT, X, Y, lat, lon, cfg, circ, jack, rot, ...
                     divArea, cxi, cyi, prof, oseen, f)
% MAKEFIGURES  All output figures.
%
%   Color choices are deliberate, not decorative:
%   * The three vorticity estimators get three DISTINCT HUES, because they are
%     different things (categorical), and each is labeled so identity never
%     rests on color alone.
%   * The eddy-frame scatter is colored by elapsed time, a MAGNITUDE, so it uses
%     a single-hue light-to-dark ramp.  A rainbow map (jet, turbo) would imply
%     categories that do not exist and is not colorblind-safe.
arguments (Input)
    t datetime
    fitTT timetable
    X double
    Y double
    lat double
    lon double
    cfg struct
    circ struct
    jack struct
    rot struct
    divArea (:,1) double
    cxi (:,1) double
    cyi (:,1) double
    prof struct
    oseen struct
    f (1,1) double
end % arguments Input

C = [0.165 0.471 0.839;      % blue   - circulation / Stokes
     0.922 0.408 0.204;      % orange - least squares
     0.106 0.686 0.478];     % green  - constellation rotation

% --- vorticity ----------------------------------------------------------
fig = figure(Name="vorticity", Position=[100 100 1000 660], Color="w");
tl = tiledlayout(fig, 4, 1, TileSpacing="compact", Padding="compact");

ax = nexttile(tl, 1, [3 1]);
hold(ax, "on"); grid(ax, "on");
lo = jack.zeta - jack.spread;  hi = jack.zeta + jack.spread;
ok = isfinite(lo) & isfinite(hi);
fill(ax, [t(ok); flipud(t(ok))], 1e3*[lo(ok); flipud(hi(ok))], C(1,:), ...
    FaceAlpha=0.16, EdgeColor="none", DisplayName="leave-one-out spread");
plot(ax, t, 2e3*rot.omega, Color=C(3,:), LineWidth=1.0, ...
    DisplayName="constellation rotation (positions only)");
plot(ax, fitTT.Time, 1e3*fitTT.zeta, Color=C(2,:), LineWidth=1.4, ...
    DisplayName="least-squares gradient");
plot(ax, t, 1e3*circ.zeta, Color=C(1,:), LineWidth=1.8, ...
    DisplayName="circulation / Stokes (primary)");
yline(ax, 0, Color=[.4 .4 .4], HandleVisibility="off");
ylabel(ax, "\zeta  (10^{-3} s^{-1})");
ylim(ax, [-4.2 1.4]);
legend(ax, Location="northeast", FontSize=8, Box="off");
zMed = median(circ.zeta, "omitnan");
title(ax, sprintf("Peleliu tip vortex: median \\zeta = %.2e s^{-1}, Rossby %.0f", ...
    zMed, zMed/f));
fontsize(ax, 9, "points");

ax2 = nexttile(tl, 4);
plot(ax2, t, circ.area/1e6, Color=[.35 .35 .35], LineWidth=1);
grid(ax2, "on");
ylabel(ax2, "cluster area (km^2)");
linkaxes([ax ax2], "x");
exportgraphics(fig, fullfile(cfg.outDir, "vorticity.png"), Resolution=200);

% --- divergence, strain, Okubo-Weiss ------------------------------------
fig = figure(Name="kinematics", Position=[100 100 1000 740], Color="w");
tl = tiledlayout(fig, 3, 1, TileSpacing="compact", Padding="compact");

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on");
plot(ax, t, circ.delta, Color=C(1,:), LineWidth=1, DisplayName="contour flux");
plot(ax, t, divArea,    Color=C(3,:), LineWidth=1, DisplayName="d(lnA)/dt");
plot(ax, fitTT.Time, fitTT.delta, Color=C(2,:), LineWidth=1, ...
    DisplayName="least squares");
yline(ax, 0, Color=[.4 .4 .4], HandleVisibility="off");
ylabel(ax, "divergence (s^{-1})");
legend(ax, Orientation="horizontal", Location="northoutside", Box="off", FontSize=8);

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on");
plot(ax, fitTT.Time, fitTT.sigN, Color=C(1,:), LineWidth=1, DisplayName="normal");
plot(ax, fitTT.Time, fitTT.sigS, Color=C(2,:), LineWidth=1, DisplayName="shear");
plot(ax, fitTT.Time, hypot(fitTT.sigN, fitTT.sigS), "k", LineWidth=1, ...
    DisplayName="total strain");
plot(ax, fitTT.Time, abs(fitTT.zeta), Color=C(3,:), LineStyle="--", ...
    LineWidth=1, DisplayName="|\zeta|");
ylabel(ax, "strain (s^{-1})");
legend(ax, Orientation="horizontal", Location="northoutside", Box="off", FontSize=8);

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on");
plot(ax, fitTT.Time, fitTT.OW, "k", LineWidth=1);
yline(ax, 0, Color="r");
ylabel(ax, "Okubo-Weiss (s^{-2})");
text(ax, 0.01, 0.12, "OW < 0: rotation dominates strain", Units="normalized", ...
    FontSize=8, Color=C(1,:));
exportgraphics(fig, fullfile(cfg.outDir, "kinematics.png"), Resolution=200);

% --- tracks, center, and a real map --------------------------------------
fig = figure(Name="tracks", Position=[100 100 1300 480], Color="w");
tl = tiledlayout(fig, 1, 3, TileSpacing="compact", Padding="compact");

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on"); axis(ax,"equal");
plot(ax, X/1e3, Y/1e3, LineWidth=0.8);
plot(ax, cxi/1e3, cyi/1e3, "k.-", MarkerSize=4, LineWidth=0.8);
xlabel(ax,"east (km)"); ylabel(ax,"north (km)"); title(ax,"Earth frame");
legend(ax, [cfg.drifters, "eddy center"], FontSize=7, Location="best", Box="off");

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on"); axis(ax,"equal");
plot(ax, X - cxi, Y - cyi, LineWidth=0.7);
plot(ax, 0, 0, "k+", MarkerSize=12, LineWidth=1.5);
xlabel(ax,"east of center (m)"); ylabel(ax,"north of center (m)");
title(ax,"Eddy frame");

% A geographic axes puts the track where it actually happened.  The basemap
% needs internet; fall back to a plain graticule if it is unavailable.
gx = geoaxes(Parent=tl);
gx.Layout.Tile = 3;
hold(gx, "on");
% geoplot takes VECTORS, not the (nEpoch x nDrifter) matrices the rest of the
% script passes around, so plot one drifter at a time.  Looping also keeps the
% per-drifter colors consistent with the Earth-frame panel beside it.
for j = 1:size(lat, 2)
    geoplot(gx, lat(:,j), lon(:,j), LineWidth=0.9);
end % for j
try
    geobasemap(gx, "satellite");
catch
    geobasemap(gx, "none");
end % try
title(gx, "Peleliu wake");
exportgraphics(fig, fullfile(cfg.outDir, "tracks_center.png"), Resolution=200);

% --- radial profile ------------------------------------------------------
fig = figure(Name="radial", Position=[100 100 1000 420], Color="w");
tl = tiledlayout(fig, 1, 2, TileSpacing="compact", Padding="compact");
rf = linspace(min(prof.r), max(prof.r), 200).';

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on");
if all(isfinite(oseen.radiusCI))
    xregion(ax, oseen.radiusCI(1), oseen.radiusCI(2), FaceColor=C(2,:), ...
        FaceAlpha=0.10, HandleVisibility="off");
end % if all
errorbar(ax, prof.r, prof.vTheta, prof.vThetaErr, "o-", Color=C(1,:), ...
    MarkerSize=4, LineWidth=1, DisplayName="v_\theta observed");
plot(ax, prof.r, prof.vR, "s--", Color=[.55 .55 .55], MarkerSize=4, ...
    LineWidth=0.8, DisplayName="v_r (should be ~0)");
plot(ax, rf, oseen.fitobj(rf), Color=C(2,:), LineWidth=1.4, ...
    DisplayName=sprintf("Lamb-Oseen R = %.0f m", oseen.radius));
xline(ax, oseen.radius, Color=C(2,:), LineStyle=":", HandleVisibility="off");
xlabel(ax,"radius from fitted center (m)"); ylabel(ax,"velocity (m s^{-1})");
legend(ax, Location="best", FontSize=8, Box="off");
title(ax, "Azimuthal velocity profile (shaded: 95% CI on R)");

ax = nexttile(tl); hold(ax,"on"); grid(ax,"on");
plot(ax, prof.r, prof.zeta, "o-", Color=C(1,:), MarkerSize=4, LineWidth=1, ...
    DisplayName="observed \zeta(r)");
plot(ax, rf, oseen.gamma/(pi*oseen.radius^2)*exp(-(rf/oseen.radius).^2), ...
    Color=C(2,:), LineWidth=1.4, DisplayName="Lamb-Oseen");
yline(ax, median(circ.zeta,"omitnan"), Color=[.4 .4 .4], LineStyle="--", ...
    DisplayName="cluster-mean \zeta");
xlabel(ax,"radius from fitted center (m)"); ylabel(ax,"\zeta (s^{-1})");
legend(ax, Location="best", FontSize=8, Box="off");
title(ax, "Local vorticity vs radius");
exportgraphics(fig, fullfile(cfg.outDir, "radial_profile.png"), Resolution=200);

% --- eddy-frame scatter, one panel per drifter, colored by elapsed time ---
fig = figure(Name="eddyFrameScatter", Position=[100 100 920 900], Color="w");
tl = tiledlayout(fig, 2, 2, TileSpacing="compact", Padding="compact");
hrs = hours(t - t(1));
xe = (X - cxi)/1e3;  ye = (Y - cyi)/1e3;
lim = ceil(prctile(abs([xe(:); ye(:)]), 99)*4)/4;
th  = linspace(0, 2*pi, 240);

ax = gobjects(1, numel(cfg.drifters));
for j = 1:numel(cfg.drifters)
    ax(j) = nexttile(tl); hold(ax(j),"on"); grid(ax(j),"on"); axis(ax(j),"equal");
    plot(ax(j), oseen.radius/1e3*cos(th), oseen.radius/1e3*sin(th), ...
        Color=[.55 .55 .55], LineStyle="--", LineWidth=1);
    plot(ax(j), 0, 0, "+", Color=[.35 .35 .35], MarkerSize=10, LineWidth=1.4);
    scatter(ax(j), xe(:,j), ye(:,j), 12, hrs, "filled");
    xlim(ax(j), [-lim lim]); ylim(ax(j), [-lim lim]);
    title(ax(j), cfg.drifters(j), FontSize=10);
    if j > 2; xlabel(ax(j), "east of eddy center (km)"); end
    if mod(j,2) == 1; ylabel(ax(j), "north of eddy center (km)"); end
end % for j
colormap(fig, singleHueRamp(256));
cb = colorbar(ax(end));
cb.Layout.Tile = "east";
cb.Label.String = "hours from common start";
title(tl, "Drifter position in the eddy reference frame");
exportgraphics(fig, fullfile(cfg.outDir, "eddy_frame_scatter.png"), Resolution=200);
end % makeFigures

function cmap = singleHueRamp(n)
% SINGLEHUERAMP  Light-to-dark blue ramp for a SEQUENTIAL variable.
%   Elapsed time is a magnitude, and magnitudes get one hue running light to
%   dark.  Rainbow maps (jet, hsv) invent perceptual boundaries the data do not
%   have and are not colorblind-safe; parula is better but still shifts hue.
%   One hue makes "later = darker" readable at a glance.
arguments (Input)
    n (1,1) double {mustBePositive, mustBeInteger} = 256
end % arguments Input
arguments (Output)
    cmap (:,3) double {mustBeInRange(cmap,0,1)}
end % arguments Output
lo = [0.827 0.890 0.969];      % #d3e3f7
hi = [0.059 0.227 0.408];      % #0f3a68
cmap = lo + (hi - lo) .* linspace(0,1,n).';
end % singleHueRamp
