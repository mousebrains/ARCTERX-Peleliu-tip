% EXPLORE_DRIFTER_PATHS  Fit an eddy center and per-drifter orbital radii to
%                        wave-buoy GPS tracks in the Peleliu tip wake.
%
% SCIENCE CONTEXT
%   ARCTERX-2023 Wake experiment, Peleliu (Palau) tip vortex shedding.  On
%   22-23 May 2023 a single shedding event was sampled by four wave buoys
%   ("mwb" drifters) whose 2 Hz GPS tracks loop around a common eddy.
%
% MODEL
%   Each drifter's Earth-frame position is written as the position of the
%   shared eddy centroid plus a radial offset from that centroid:
%
%       \vec{x_{earth}} = \vec{x_{center}}
%                       + radius_{drifter} * \vec{\theta_{drifter}}
%
%   The orbital direction \theta is not fitted; it is taken from the
%   drifter's own smoothed velocity, assuming the drifter moves tangentially
%   around the eddy so the center lies +/-90 degrees off the velocity vector.
%   Only the centroid (x,y) and one radius per drifter are free parameters.
%   Because all four drifters orbit the SAME center, they are fitted jointly
%   in each time window: 2 + nDrifters parameters, 4 tracks of constraint.
%
% INPUT
%   <dataDir>/mwb*_gps_timeseries.nc  -- CF-1.6 trajectory files from
%   Scripps/CORDC.  Variables are stored as (time, sampling_dwell) blocks of
%   2048 samples; loadData() flattens them to a single time series.
%     time    milliseconds since 1970-01-01 (UTC)
%     u,v,w   m/s, eastward/northward/upward GPS velocity
%     lat,lon degrees north / east (WGS-84)
%     sog     m/s over ground,  cog degrees,  numsats, pdop  (diagnostics)
%
% METHOD
%   1. Load each drifter, de-duplicate and sort by time.
%   2. Retime all drifters onto one common grid (median sample interval),
%      leaving gaps as missing.
%   3. Convert lat/lon to a local east/north metric frame about the mean
%      position, then low-pass both the velocity and the position with a
%      60 s lowess window to suppress orbital wave motion and GPS noise.
%   4. In overlapping time windows, fit the shared eddy center and the four
%      radii with fminsearch (see myCentroid below).
%   5. Interpolate the fitted center back onto the full time base, express
%      each drifter in the eddy frame, and reconstruct the modeled track.
%
% OUTPUT  (written to outDir)
%   fits.mat        tbl (per-drifter timetables) and centroids (fit table)
%   tracks.png      Earth-frame drifter tracks + fitted eddy center path
%   eddy_frame.png  drifter positions relative to the fitted eddy center
%   radius.png      fitted orbital radius of each drifter versus time
%   fit.png         Earth-frame tracks with the reconstructed model tracks
%
% REQUIREMENTS
%   MATLAB R2023b or later (smoothdata2), Mapping Toolbox (distance,
%   wgs84Ellipsoid).  Parallel Computing Toolbox is optional -- without it
%   the parfor loops run serially.
%
% CAVEATS  (see also the NOTES section at the end of this file)
%   * The fit is unconstrained, so a drifter's radius may come out negative;
%     that is the same degeneracy as flipping which side of the velocity
%     vector the center sits on.  The radius figure plots abs(radius).
%   * fminsearch stops on whichever of MaxIter / MaxFunEvals is reached
%     first.  MaxFunEvals is left at its default (200*nPars), so it, not the
%     larger MaxIter set below, is the binding limit.  Windows that stop on
%     a limit return exitflag == 0 and are dropped from the centroid series.
%
% July-2026, Pat Welch, pat@mousebrains.com

%% ------------------------------------------------------------------------
%  Configuration
%  ------------------------------------------------------------------------

% Drifters deployed on the date of interest
drifters = append("mwb", ["458d02", "788d01", "790d01", "793d02"]);

myPath  = fileparts(mfilename("fullpath")); % Script's path
dataDir = myPath;                           % Where the NetCDF files live
outDir  = myPath;                           % Where figures and fits.mat go

smoothWidth = seconds(60);  % Low-pass width for velocity and position
dtStep      = minutes(5);   % Time between the start of successive fit windows
dtWindow    = 4 * dtStep;   % Width of each fit window (so windows overlap 4x)
radius0     = 1000;         % Initial guess at each drifter's orbit radius (m)

%% ------------------------------------------------------------------------
%  Load the drifters
%  ------------------------------------------------------------------------

% Get filenames of NetCDF files to load
files = struct2table(dir(fullfile(dataDir, "*.nc")));
files.title = extractBefore(string(files.name), "_");
files = files(ismember(files.title, drifters),:);
files.fn = string(fullfile(files.folder, files.name));

if isempty(files)
    error("No files matching %s found in %s", strjoin(drifters, ", "), dataDir);
end % if isempty

% Load the NetCDF files into tbl, one cell/file
nDrifters = size(files,1);
tbl = cell(nDrifters, 1);

t0 = NaT; % Earliest observation across all drifters
t1 = NaT; % Latest observation across all drifters
dtDrifter = nan(nDrifters,1); % Median sample interval of each drifter (seconds)

for index = 1:nDrifters
    tbl{index} = loadData(files.fn{index});
    t0 = min(tbl{index}.time(1), t0);
    t1 = max(tbl{index}.time(end), t1);
    dtDrifter(index) = seconds(median(diff(tbl{index}.time)));
end % for index

dt = seconds(median(dtDrifter, "omitmissing")); % Common sample interval

%% ------------------------------------------------------------------------
%  Put every drifter onto one common time base
%  ------------------------------------------------------------------------

% Time to project everything onto.  retime defaults to "fillwithmissing", so
% gaps stay missing rather than being interpolated across.
tJoint = t0:dt:t1;
tbl = cellfun(@(x) retime(x, tJoint), tbl, UniformOutput=false);

% Center point for distances
lat0 = mean(cellfun(@(x) mean(x.lat, "omitmissing"), tbl));
lon0 = mean(cellfun(@(x) mean(x.lon, "omitmissing"), tbl));

% Smooth over wave motion and GPS noise.  The drifters orbit on the waves at
% ~5-15 s, so a 60 s window removes the wave signal while leaving the eddy.
windowWidth = round(smoothWidth / dt); % Samples, 120 at 2 Hz

fprintf("%d drifters, %s to %s, dt = %s, smoothing window = %d samples\n", ...
    nDrifters, t0, t1, dt, windowWidth);

error("Gotme");
%% ------------------------------------------------------------------------
%  Smooth, and build the local east/north frame
%  ------------------------------------------------------------------------

tStart = tic();
parfor index = 1:nDrifters
    fprintf("Smoothing %d\n", index);
    a = tbl{index};

    % Smooth velocity to get the orbital direction.  smoothdata2 treats
    % [u v] as a 2-D field, but the second dimension is only 2 wide, so away
    % from the ends of the record this matches per-column smoothing to
    % ~1e-6 m/s; only the first/last window-width of samples differ.
    su = smoothdata2([a.u, a.v], "lowess", windowWidth, "includemissing");
    a.su = su(:,1);
    a.sv = su(:,2);

    % Angle of the velocity vector in the cartesian plane, degrees CCW from east
    a.theta_velocity = atan2d(a.sv, a.su);

    % Geodesic range and azimuth from the reference point, converted to a
    % local east/north frame.  az is degrees clockwise from north, so
    % east = dist*cosd(90-az) and north = dist*sind(90-az).
    [dist, az] = distance(lat0, lon0, a.lat, a.lon, wgs84Ellipsoid("meter"));
    a.x = dist .* cosd(90-az); % Meters east of (lat0, lon0)
    a.y = dist .* sind(90-az); % Meters north of (lat0, lon0)

    % Same low-pass applied to position
    sx = smoothdata2([a.x, a.y], "lowess", windowWidth, "includemissing");
    a.sx = sx(:,1);
    a.sy = sx(:,2);

    tbl{index} = a;
end % parfor index

fprintf("%.2f seconds for smoothing\n", toc(tStart));

%% ------------------------------------------------------------------------
%  Fit the shared eddy center and per-drifter radii, window by window
%  ------------------------------------------------------------------------
%  Since all the drifters share the same eddy center, we can fit them
%  together with different radii.  Windows overlap (dtWindow = 4*dtStep) so
%  the centroid series is smoother than the individual fits.

stimeWindows = t0:dtStep:t1; % Starting time of each window
centroids = cell(size(stimeWindows));

nPars = 2 + nDrifters; % xCenter, yCenter, one radius per drifter
nTimeWindows = numel(stimeWindows);

tStart = tic();
parfor sindex = 1:nTimeWindows
    startTime = tic();
    stime = stimeWindows(sindex);
    etime = stime + dtWindow;

    % Pull the rows of every drifter that fall inside this window
    b = cell(nDrifters,1); % For pruned dataset
    for index = 1:nDrifters
        a = tbl{index}(:,["theta_velocity", "sx", "sy"]);
        a = rmmissing(a(a.time >= stime & a.time < etime,:));
        a.index = repmat(index, size(a.time));
        b{index} = a;
    end % for index
    b = vertcat(b{:});

    if isempty(b) % Nothing to fit, e.g. a gap in every drifter
        fprintf("%s to %s has no data\n", stime, etime);
        continue;
    end % if isempty

    pars0 = [mean(b.sx), mean(b.sy), repmat(radius0, 1, nDrifters)];
    [pars, ~, exitflag, output] = fminsearch(@(x) myCentroid(x, b), pars0, ...
        optimset(Display="notify", MaxFunEvals=4000*nPars, MaxIter=4000*nPars));

    if exitflag % Success, i.e. converged rather than hitting a limit
        acent = table();
        acent.time = median(b.time); % Middle of the data actually used
        acent.xCenter = pars(1);
        acent.yCenter = pars(2);
        for index = 1:nDrifters
            acent.(sprintf("drifter_%d", index)) = pars(index + 2);
        end % for index
        centroids{sindex} = acent;
    else
        fprintf("%s failed\n", stime);
        disp(output)
    end % if exitflag

    fprintf("%.2f seconds for %s to %s\n", toc(startTime), stime, etime);
end % parfor sindex

% parfor fills centroids by window index, so the surviving rows are already
% in time order, which interp1 below relies on.
centroids = centroids(~cellfun(@isempty, centroids)); % Drop empty cells
centroids = vertcat(centroids{:});

fprintf("%.2f seconds for %d of %d fits\n", ...
    toc(tStart), size(centroids,1), numel(stimeWindows));

%% ------------------------------------------------------------------------
%  Use the projected center to recalculate radius for each drifter in
%  the time window. This takes into account the moving eddy during the
%  window duration.
%  ------------------------------------------------------------------------

radii = nan(nTimes, nDrifters); % A radius for each time window and drifter
for index = 1:nDrifters
    startTime = tic();
    a = tbl{index}; %  drifter's information
    name = sprintf("drifter_%d", index); % centroid column name
    % Center at each instant in time, C2
    a.xCenter = interp1(centroid.time, centroid.xCenter, a.time, "cubic");
    a.yCenter = interp1(centroid.time, centroid.yCenter, a.time, "cubic"); 
    a.dx = a.sx - a.xCenter;
    a.dy = a.sy - a.yCenter;
    a = a(~ismissing(a.dx) & ~isnan(a.dy),:)
end % for index

%% ------------------------------------------------------------------------
%  Project each drifter into the eddy frame and rebuild the modeled track
%  ------------------------------------------------------------------------

tStart = tic();
parfor index = 1:nDrifters
    a = tbl{index};

    % Calculate drifter reference frame
    a.xCenter = interp1(centroids.time, centroids.xCenter, a.time, "linear");
    a.yCenter = interp1(centroids.time, centroids.yCenter, a.time, "linear");
    a.xEddy = a.sx - a.xCenter; % Meters east of the eddy center
    a.yEddy = a.sy - a.yCenter; % Meters north of the eddy center

    % Fitted track: center plus the fitted radius along the direction
    % perpendicular to the drifter's velocity
    name = sprintf("drifter_%d", index);
    radius = interp1(centroids.time, centroids.(name), a.time, "linear");
    a.xFit = a.xCenter + radius .* cosd(a.theta_velocity + 90);
    a.yFit = a.yCenter + radius .* sind(a.theta_velocity + 90);

    tbl{index} = a;
end % parfor index
fprintf("%.2f seconds for projection\n", toc(tStart));

save(fullfile(outDir, "fits"), "tbl", "centroids");

%% ------------------------------------------------------------------------
%  Figures
%  ------------------------------------------------------------------------

colors = jet(nDrifters);
tTitle = sprintf("%s to %s UTC", t0, t1);

% Figure 1: Earth-frame tracks, raw and smoothed, with the eddy center path
figure(1);
clf;
h = gobjects(nDrifters,1);
for index = 1:nDrifters
    a = tbl{index};
    h(index) = scatter(a.x/1000, a.y/1000, 10, colors(index,:));
    hold on;
    plot(a.sx/1000, a.sy/1000, "-", Color=colors(index,:));
end % for index
hc = plot(centroids.xCenter/1000, centroids.yCenter/1000, "ko-", MarkerSize=3);
hold off;
axis tight;
daspect([1,1,1]);
grid on;
xlabel(sprintf("Eastward from %.6f degrees (km)", lon0));
ylabel(sprintf("Northward from %.6f degrees (km)", lat0));
legend([h; hc], [files.title; "eddy center"], Location="best");
title(tTitle);
print(fullfile(outDir, "tracks.png"), "-dpng");

% Figure 2: positions relative to the fitted eddy center
figure(2);
clf;
t = tiledlayout(max(1, floor(sqrt(nDrifters))), max(1, ceil(sqrt(nDrifters))), ...
    TileSpacing="tight", Padding="tight");
h = gobjects(prod(t.GridSize), 1);
for index = 1:nDrifters
    a = tbl{index};
    h(index) = nexttile(index);
    scatter(a.xEddy/1000, a.yEddy/1000, 10, hours(a.time - t0));
    grid on;
    daspect([1,1,1]);
    axis tight;
    legend(files.title(index), Location="best");
    if rem(index,2) == 0 % even
        ylabel("Northward from eddy center (km)");
    else % odd
        xlabel("Eastward from eddy center (km)");
    end % if rem
end % for index

cb = colorbar();
cb.Layout.Tile = "south";
cb.Label.String = sprintf("Hours since %s", t0);

sgtitle(tTitle);
print(fullfile(outDir, "eddy_frame.png"), "-dpng");

% Figure 3: fitted orbital radius versus time.  abs() because the fit is
% unconstrained and the sign only says which side of the velocity the
% center was placed on.
figure(3);
clf;
h = gobjects(nDrifters,1);
for index = 1:nDrifters
    name = sprintf("drifter_%d", index);
    h(index) = plot(centroids.time, abs(centroids.(name)) / 1000, "-", ...
        Color=colors(index,:));
    hold on;
end % for index
hold off
axis tight;
grid on;
xlabel("Time (UTC)");
ylabel("Radius (km)");
legend(h, files.title, Location="best");
title(tTitle);
print(fullfile(outDir, "radius.png"), "-dpng");

% Figure 4: as figure 1, plus the track reconstructed from the fit
figure(4);
clf;
h = gobjects(nDrifters,1);
for index = 1:nDrifters
    a = tbl{index};
    h(index) = scatter(a.x/1000, a.y/1000, 10, colors(index,:));
    hold on;
    plot(a.sx/1000, a.sy/1000, "-", Color=colors(index,:));
    plot(a.xFit/1000, a.yFit/1000, "--", Color=colors(index,:));
end % for index
hc = plot(centroids.xCenter/1000, centroids.yCenter/1000, "ko-", MarkerSize=3);
hold off;
axis tight;
daspect([1,1,1]);
grid on;
xlabel(sprintf("Eastward from %.6f degrees (km)", lon0));
ylabel(sprintf("Northward from %.6f degrees (km)", lat0));
legend([h; hc], [files.title; "eddy center"], Location="best");
title(tTitle);
print(fullfile(outDir, "fit.png"), "-dpng");

%% ------------------------------------------------------------------------
%  NOTES / future work
%  ------------------------------------------------------------------------
%
% We can write the position of a drifter at time t as
% the centroid of the eddy plus a position on
% an evolving circle about the centroid
%
% x_drifter(t) = x_centroid(t) + x'_drifter(t)
%   where these are vectors and
%   x' = r_drifter(t) * cosine(theta_drifter(t))
%   y' = r_drifter(t) * sine(theta_drifter(t))
%
% We have physical constraints about how much
% the centroid and drifters can move in a given timestep
% due to continuity.
%
% dx_centroid(t) = x_centroid(t+epsilon) - x_centroid(t-epsilon)
% dr_drifter(t) = r_drifter(t+epsilon) - r_drifter(t-epsilon)
% dtheta_drifter(t) = theta_drifter(t+epsilon) - theta_drifter(t-epsilon)
%
% Theta has a directional constraint too, dtheta_drifter should be
% consistently one sign.
%
% Windows are currently fitted independently, so none of the above continuity
% constraints are imposed and successive centroids can jump.  A joint fit
% over all windows, or a penalty on d(centroid)/dt, would tie them together.

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function dist = myCentroid(pars, b)
% MYCENTROID  Cost function for the joint eddy-center / radii fit.
%
%   Given a trial centroid and one radius per drifter, step the drifter
%   towards where the center should be -- a distance `radius` along the
%   direction perpendicular to the drifter's velocity -- and accumulate how
%   far that estimate lands from the trial centroid.  The cost is the sum of
%   Euclidean distances (an L1 norm over the residual lengths), which is
%   less sensitive to outliers than a sum of squares.
%
%   The perpendicular direction is ambiguous: the center could be 90 degrees
%   to either side of the velocity.  Both signs are evaluated per drifter and
%   the closer one is kept, so the fit does not have to guess the rotation
%   sense.  Note this is decided independently for each drifter and each
%   window, and it is redundant with the sign of the fitted radius, so the
%   returned radii can be negative.
%
%   Inputs
%     pars  [xCentroid; yCentroid; radius_1; ... radius_n], meters.
%           fminsearch passes back a vector with the shape of the initial
%           guess (a row here); the arguments block reshapes it to a column.
%     b     Timetable of the window's samples, with sx, sy (meters),
%           theta_velocity (degrees) and index (which drifter the row is).
%
%   Output
%     dist  Total distance between the observed and estimated center (m).

arguments (Input)
    pars (:,1) double
    b timetable
end % arguments Input
arguments (Output)
    dist (1,1) double % distance from observed drifter position and estimated position
end % arguments Output

xCentroid = pars(1);
yCentroid = pars(2);
radii = pars(3:end);

dist = 0;
for index = 1:numel(radii)
    radius = radii(index);
    bb = b(b.index == index,:);
    if isempty(bb), continue; end
    theta_perp = bb.theta_velocity + 90; % Perpendicular to velocity direction
    rx = radius * cosd(theta_perp);
    ry = radius * sind(theta_perp);

    dxp = sum(sqrt((xCentroid - (bb.sx + rx)).^2 + (yCentroid - (bb.sy + ry)).^2));
    dxm = sum(sqrt((xCentroid - (bb.sx - rx)).^2 + (yCentroid - (bb.sy - ry)).^2));
    dx = min(dxp, dxm); % Pick sign closest to zero
    dist = dist + dx;
end % for index
end % myCentroid

function tbl = loadData(fn)
% LOADDATA  Read one CORDC wave-buoy GPS NetCDF file into a timetable.
%
%   The file stores each variable as (time, sampling_dwell), i.e. blocks of
%   2048 samples per wave-processing period, so every variable is flattened
%   to a single column.  Duplicate timestamps are dropped and the record is
%   returned in time order (unique sorts as a side effect).  Times in the
%   file are milliseconds since 1970-01-01 UTC.

arguments (Input)
    fn string {mustBeFile}
end % arguments Input
arguments (Output)
    tbl timetable
end % arguments Output

names = ["time", "u", "v", "w", "lat", "lon", "sog", "cog", "numsats", "pdop"];

fprintf("Reading %s\n", fn);

tbl = table();
for name = names
    val = ncread(fn, name);
    tbl.(name) = val(:);
end % for name

[~, ix] = unique(tbl.time); % Also sorts into time order
tbl = tbl(ix,:);
tbl.spd = sqrt(tbl.u.^2 + tbl.v.^2); % Horizontal speed, diagnostic only

tbl.time = datetime(tbl.time / 1000, ConvertFrom="posixtime");
tbl = table2timetable(tbl, RowTimes="time");
end % loadData
