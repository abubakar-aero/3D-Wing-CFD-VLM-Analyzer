%% ========================================================================
%  3D VORTEX LATTICE METHOD (VLM) - INTERACTIVE WING AERODYNAMIC ANALYSIS
%  ------------------------------------------------------------------------
%  A self-contained MATLAB tool for low-order potential-flow analysis of
%  custom wing planforms using discrete horseshoe vortices.
%
%  MODELING ASSUMPTIONS (standard for classical VLM):
%   - Thin, zero-thickness lifting surface (flat panels, no camber/twist
%     unless added by the user).
%   - Inviscid, incompressible, irrotational flow (potential flow). No
%     viscous drag, no compressibility, no stall modeling.
%   - Steady-state solution only.
%   - Trailing vortices are assumed to leave the trailing edge parallel to
%     the freestream direction and extend to a numerically "infinite"
%     distance downstream.
%
%  Author: Principal Aerospace Engineer / Senior MATLAB Developer (demo)
% =========================================================================
clear; clc; close all;

%% ================= 1. INTERACTIVE POP-UP CONFIGURATION ==================
% Fixed (non-dialog) defaults - dihedral, air density and mesh density are
% not exposed in the pop-ups per spec, but can be edited here directly.
config.dihedral_deg = 3;       % Dihedral angle [deg]
config.rho          = 1.225;   % Air density, rho [kg/m^3]
config.n_span        = 15;      % Spanwise panel divisions PER SIDE
config.n_chord        = 8;      % Chordwise panel divisions

validWingTypes = {'RECTANGULAR','SWEPT','DELTA','FACETED_STEALTH','BLENDED_WING_BODY'};

config = get_config_via_dialogs(config, validWingTypes);

fprintf('\n--- Running VLM analysis: %s wing, %d x %d panels ---\n', ...
    config.wing_type, 2*config.n_span, config.n_chord);

%% ================= 2. FREESTREAM & REFERENCE QUANTITIES ================
alpha_rad = deg2rad(config.alpha_deg);
% Body axes convention: x = forward/aft (chordwise), y = spanwise,
% z = up. Freestream vector tilted by angle of attack in the x-z plane.
Vinf_vec = config.V_inf * [cos(alpha_rad), 0, sin(alpha_rad)];
q_inf    = 0.5 * config.rho * config.V_inf^2;   % Dynamic pressure

%% ========================= 3. GEOMETRY / MESH ===========================
tic;
panels = build_wing_panels(config);
N = numel(panels);
fprintf('Mesh generated: N = %d panels (%.2f s)\n', N, toc);

%% ==================== 4. INFLUENCE MATRIX & RHS =========================
tic;
[Amat, RHS] = assemble_system(panels, Vinf_vec, config);
fprintf('Influence matrix assembled: %d x %d (%.2f s)\n', N, N, toc);

%% ============================ 5. SOLVE ==================================
Gamma = Amat \ RHS;    % Circulation strength of each horseshoe vortex [m^2/s]

%% ======================= 6. AERODYNAMIC FORCES ==========================
results = compute_aero_forces(panels, Gamma, config, Vinf_vec, q_inf);

fprintf('\n=== RESULTS ===\n');
fprintf('S_ref = %.3f m^2\n', results.S_ref);
fprintf('Lift  = %.1f N\n', results.L_total);
fprintf('C_L   = %.4f\n', results.CL);
fprintf('C_Di  = %.5f\n', results.CDi);
fprintf('L/D   = %.2f\n', results.LD);

%% ========================= 7. VISUALIZATION =============================
plot_dashboard(panels, Gamma, results, config);


%% ========================================================================
%                          LOCAL FUNCTIONS
% =========================================================================

function config = get_config_via_dialogs(config, validWingTypes)
% Pop-up interactive setup: a menu() selector for wing planform, followed
% by an inputdlg() multi-field box for the flight/geometry parameters.
% Cancelling either dialog aborts the script with a clear error rather
% than silently proceeding with a partially-configured analysis.

    %% --- Wing selection pop-up menu -------------------------------------
    choice = menu('Select Wing Planform', validWingTypes{:});
    if choice == 0
        error('VLM:UserCancelled', ...
            'Wing selection was cancelled/closed. Script aborted.');
    end
    config.wing_type = validWingTypes{choice};

    %% --- Flight/geometry parameter input dialog -------------------------
    prompt_labels = { ...
        'Wingspan  b  (m):', ...
        'Root chord  c_root  (m):', ...
        'Tip chord  c_tip  (m):', ...
        'LE sweep angle  \Lambda  (deg):', ...
        'Angle of attack  \alpha  (deg):', ...
        'Freestream velocity  V_inf  (m/s):'};
    dlg_title   = sprintf('VLM Parameters - %s Wing', config.wing_type);
    num_lines   = 1;
    defaults    = {'10.0', '2.5', '0.8', '30', '6', '50'};

    answer = inputdlg(prompt_labels, dlg_title, num_lines, defaults);

    if isempty(answer)
        error('VLM:UserCancelled', ...
            'Parameter input was cancelled/closed. Script aborted.');
    end

    values = str2double(answer);
    if any(isnan(values))
        badIdx = find(isnan(values));
        error('VLM:InvalidInput', ...
            'Non-numeric entry for: %s. Please re-run and enter numeric values.', ...
            strjoin(prompt_labels(badIdx), ', '));
    end

    config.span            = values(1);
    config.root_chord       = values(2);
    config.tip_chord        = values(3);
    config.sweep_angle_deg  = values(4);
    config.alpha_deg        = values(5);
    config.V_inf             = values(6);

    % --- Sanity checks on physically meaningful ranges -------------------
    if config.span <= 0 || config.root_chord <= 0 || config.tip_chord < 0
        error('VLM:InvalidInput', 'Span and chords must be positive values.');
    end
    if config.V_inf <= 0
        error('VLM:InvalidInput', 'Freestream velocity must be positive.');
    end
end

% ------------------------------------------------------------------------
%  GEOMETRY / PANEL MESH GENERATOR
% ------------------------------------------------------------------------
function [chord_fn, LE_fn] = get_planform_functions(config)
% Returns vectorized function handles chord(y) and LE_x(y), y in [-b/2,b/2]

    b     = config.span;
    croot = config.root_chord;
    ctip  = config.tip_chord;
    sweep = deg2rad(config.sweep_angle_deg);
    halfb = b/2;

    switch config.wing_type
        case 'RECTANGULAR'
            % Constant chord across span: c(y) = c_root
            chord_fn = @(y) croot * ones(size(y));
            LE_fn    = @(y) zeros(size(y));

        case 'SWEPT'
            % Linear taper + straight-line LE sweep
            chord_fn = @(y) croot - (croot - ctip) .* (abs(y) / halfb);
            LE_fn    = @(y) abs(y) * tan(sweep);

        case 'DELTA'
            % Swept LE, straight (unswept) TE -> triangular planform that
            % tapers to a sharp point near the tip.
            LE_fn    = @(y) abs(y) * tan(sweep);
            chord_fn = @(y) max(croot - LE_fn(y), 0.02*croot);

        case 'FACETED_STEALTH'
            % Kinked ("double-sweep") leading edge: two straight facet
            % segments per half-span, producing a faceted trapezoidal
            % planform reminiscent of low-observable airframe outlines.
            kink        = 0.55 * halfb;
            sweep_outer = deg2rad(max(config.sweep_angle_deg - 20, 5));
            LE_fn    = @(y) faceted_LE(y, kink, sweep, sweep_outer);
            chord_fn = @(y) croot - (croot - ctip) .* (abs(y) / halfb);

        case 'BLENDED_WING_BODY'
            % Wide, constant-chord center body smoothly (cosine) blended
            % into a swept outer wing panel.
            blend    = 0.30 * halfb;
            LE_fn    = @(y) bwb_LE(y, blend, sweep);
            chord_fn = @(y) bwb_chord(y, blend, halfb, croot, ctip);
    end
end

function LE = faceted_LE(y, kink, sweepIn, sweepOut)
    ay = abs(y);
    LE = zeros(size(ay));
    inner = ay <= kink;
    LE(inner)  = ay(inner) * tan(sweepIn);
    LE_kink    = kink * tan(sweepIn);
    LE(~inner) = LE_kink + (ay(~inner) - kink) * tan(sweepOut);
end

function LE = bwb_LE(y, blend, sweep)
    ay = abs(y);
    LE = zeros(size(ay));
    center = ay <= blend;
    LE(center)  = 0;                                   % unswept center body
    LE(~center) = (ay(~center) - blend) * tan(sweep);   % swept outboard panel
end

function c = bwb_chord(y, blend, halfb, croot, ctip)
    ay = abs(y);
    c = zeros(size(ay));
    center = ay <= blend;
    c(center) = croot;                                  % wide center body
    s = (ay(~center) - blend) / (halfb - blend);         % 0..1
    s = 0.5 - 0.5*cos(pi*s);                              % cosine smoothing
    c(~center) = croot + s .* (ctip - croot);
end

function panels = build_wing_panels(config)
% Builds the full (both-sides) panel mesh: N = 2*n_span*n_chord panels.
% Each panel stores its 4 corner vertices, quarter-chord bound-vortex
% endpoints, three-quarter-chord collocation point, outward unit normal,
% panel area, and spanwise width dy.

    b        = config.span;
    n_span   = config.n_span;
    n_chord  = config.n_chord;
    dihedral = deg2rad(config.dihedral_deg);

    [chord_fn, LE_fn] = get_planform_functions(config);

    y_edges = linspace(-b/2, b/2, 2*n_span + 1);
    z_edges = abs(y_edges) * tan(dihedral);
    c_edges = chord_fn(y_edges);
    le_edges = LE_fn(y_edges);

    nPanelsSpan = 2*n_span;
    N = nPanelsSpan * n_chord;

    panels(N) = struct('P1',[],'P2',[],'P3',[],'P4',[], ...
                        'boundA',[],'boundB',[],'colloc',[], ...
                        'normal',[],'area',[],'dy',[]);

    idx = 0;
    for s = 1:nPanelsSpan
        y1 = y_edges(s);   y2 = y_edges(s+1);
        z1 = z_edges(s);   z2 = z_edges(s+1);
        le1 = le_edges(s); le2 = le_edges(s+1);
        c1 = c_edges(s);   c2 = c_edges(s+1);

        for c = 1:n_chord
            f0 = (c-1)/n_chord;
            f1 = c/n_chord;

            P1 = [le1 + f0*c1, y1, z1];   % inboard,  LE-side corner
            P4 = [le1 + f1*c1, y1, z1];   % inboard,  TE-side corner
            P2 = [le2 + f0*c2, y2, z2];   % outboard, LE-side corner
            P3 = [le2 + f1*c2, y2, z2];   % outboard, TE-side corner

            % Bound vortex placed at this panel's LOCAL quarter-chord line
            A_pt = P1 + 0.25*(P4 - P1);
            B_pt = P2 + 0.25*(P3 - P2);

            % Collocation point at this panel's LOCAL three-quarter chord
            C1 = P1 + 0.75*(P4 - P1);
            C2 = P2 + 0.75*(P3 - P2);
            colloc = 0.5*(C1 + C2);

            % Outward normal & area via quadrilateral diagonals
            diag1 = P3 - P1;
            diag2 = P2 - P4;
            crossD = cross(diag1, diag2);
            area = 0.5 * norm(crossD);
            nvec = crossD / norm(crossD);
            if nvec(3) < 0
                nvec = -nvec;   % enforce "upward" (+z-leaning) convention
            end

            idx = idx + 1;
            panels(idx).P1 = P1;  panels(idx).P2 = P2;
            panels(idx).P3 = P3;  panels(idx).P4 = P4;
            panels(idx).boundA = A_pt;
            panels(idx).boundB = B_pt;
            panels(idx).colloc = colloc;
            panels(idx).normal = nvec;
            panels(idx).area   = area;
            panels(idx).dy     = (y2 - y1);
        end
    end
end

% ------------------------------------------------------------------------
%  VORTEX LATTICE MATHEMATICAL ENGINE
% ------------------------------------------------------------------------
function M = collect_field(panels, field)
% Stacks a 1x3 struct-array field into an Nx3 matrix.
    M = reshape([panels.(field)], 3, numel(panels))';
end

function V = biot_savart_vec(r, V1, V2)
% Vectorized Biot-Savart law: induced velocity at point r (1x3) due to a
% set of M straight vortex filaments of unit circulation strength, each
% running from V1(k,:) to V2(k,:).  Returns Mx3 induced velocities.
%
%   V = (1/4*pi) * (r1 x r2)/|r1 x r2|^2 * r0.(r1/|r1| - r2/|r2|)
%
% Singularity protection (eps = 1e-6) guards against division by zero
% when the evaluation point lies on or very near the vortex filament.
    eps_sing = 1e-6;

    R1 = r - V1;              % Mx3
    R2 = r - V2;               % Mx3
    R0 = V2 - V1;               % Mx3

    C   = cross(R1, R2, 2);      % Mx3
    Csq = sum(C.^2, 2);           % Mx1

    n1 = sqrt(sum(R1.^2, 2));
    n2 = sqrt(sum(R2.^2, 2));

    bad = (Csq < eps_sing) | (n1 < eps_sing) | (n2 < eps_sing);
    Csq(bad) = 1; n1(bad) = 1; n2(bad) = 1;   % placeholder to avoid /0

    dotterm = sum(R0 .* (R1./n1 - R2./n2), 2);   % Mx1
    K = dotterm ./ (4*pi*Csq);                     % Mx1
    V = K .* C;                                     % Mx3

    V(bad, :) = 0;   % zero out singular contributions
end

function [Amat, RHS] = assemble_system(panels, Vinf_vec, config)
% Builds the N x N influence coefficient matrix A (A_ij = induced normal
% velocity at collocation point i from a unit-strength horseshoe on panel
% j) and the RHS vector B_i = -Vinf . n_hat_i (flow-tangency condition).
%
% Each horseshoe vortex is modeled as three straight segments of unit
% circulation: an upstream trailing leg (from a numerically "infinite"
% downstream point into bound point A), the bound segment (A -> B) at the
% panel's local quarter-chord, and a downstream trailing leg (B -> a
% numerically "infinite" downstream point) - all parallel to freestream.

    N = numel(panels);
    Vinf_hat = Vinf_vec / norm(Vinf_vec);
    Lref = 1e4 * config.span;   % "practically infinite" trailing-leg length

    boundA  = collect_field(panels, 'boundA');   % Nx3
    boundB  = collect_field(panels, 'boundB');   % Nx3
    normals = collect_field(panels, 'normal');   % Nx3
    colloc  = collect_field(panels, 'colloc');   % Nx3

    A_inf = boundA + Lref*Vinf_hat;
    B_inf = boundB + Lref*Vinf_hat;

    Amat = zeros(N, N);
    RHS  = zeros(N, 1);

    for i = 1:N
        P = colloc(i, :);
        V1 = biot_savart_vec(P, A_inf, boundA);   % upstream trailing legs, all j
        V2 = biot_savart_vec(P, boundA, boundB);  % bound vortices, all j
        V3 = biot_savart_vec(P, boundB, B_inf);   % downstream trailing legs, all j
        Vij = V1 + V2 + V3;                        % Nx3: induced by each panel j

        Amat(i, :) = Vij * normals(i, :)';
        RHS(i) = -dot(Vinf_vec, normals(i, :));    % B_i = -Vinf . n_hat_i
    end
end

% ------------------------------------------------------------------------
%  AERODYNAMIC FORCE POST-PROCESSING
% ------------------------------------------------------------------------
function results = compute_aero_forces(panels, Gamma, config, Vinf_vec, q_inf)
    N    = numel(panels);
    Vinf = norm(Vinf_vec);
    Vinf_hat = Vinf_vec / Vinf;
    Lref = 1e4 * config.span;
    rho  = config.rho;

    boundA = collect_field(panels, 'boundA');
    boundB = collect_field(panels, 'boundB');
    A_inf  = boundA + Lref*Vinf_hat;
    B_inf  = boundB + Lref*Vinf_hat;

    dy = [panels.dy]';
    S  = [panels.area]';
    y_pos = 0.5*(boundA(:,2) + boundB(:,2));

    % Local pressure coefficient jump: dCp_i = 2*Gamma_i*dy_i / (Vinf*S_i)
    dCp = 2*Gamma.*dy ./ (Vinf*S + eps);

    % Sectional lift (Kutta-Joukowski) per horseshoe: dL_i = rho*Vinf*Gamma_i*dy_i
    dL = rho * Vinf * Gamma .* dy;

    % Induced drag via near-field downwash: evaluate the velocity induced
    % by every horseshoe's TRAILING legs only (bound segment excluded) at
    % each panel's own bound-vortex midpoint, then apply a rotated
    % Kutta-Joukowski relation: dDi_i = -rho * w_i * Gamma_i * dy_i
    Pmid = 0.5*(boundA + boundB);
    w = zeros(N, 1);
    for i = 1:N
        Vtrail = biot_savart_vec(Pmid(i,:), A_inf, boundA) + ...
                 biot_savart_vec(Pmid(i,:), boundB, B_inf);
        w(i) = sum(Gamma .* Vtrail(:,3));
    end
    dDi = -rho * w .* Gamma .* dy;

    S_ref    = sum(S);
    L_total  = sum(dL);
    Di_total = sum(dDi);

    CL  = L_total  / (q_inf*S_ref);
    CDi = Di_total / (q_inf*S_ref);
    LD  = CL / max(abs(CDi), 1e-9);

    % Spanwise lift distribution (aggregate chordwise strips per station)
    [y_unique, ~, ic] = unique(round(y_pos, 6));
    lift_per_station = accumarray(ic, dL);
    dy_per_station    = accumarray(ic, dy) ./ accumarray(ic, ones(N,1));
    lift_per_unit_span = lift_per_station ./ dy_per_station;

    results.dCp                = dCp;
    results.dL                 = dL;
    results.L_total             = L_total;
    results.Di_total            = Di_total;
    results.CL                  = CL;
    results.CDi                 = CDi;
    results.LD                  = LD;
    results.S_ref               = S_ref;
    results.q_inf               = q_inf;
    results.y_span              = y_unique;
    results.lift_per_unit_span  = lift_per_unit_span;
end

% ------------------------------------------------------------------------
%  FLOW-FIELD SAMPLING (for streamlines)
% ------------------------------------------------------------------------
function [X, Y, Z, U, V, W] = build_flow_grid(panels, Gamma, config)
% Samples the total velocity field (freestream + induced) on a structured
% 3D grid around the wing, for use with MATLAB's streamline() function.

    b     = config.span;
    croot = config.root_chord;
    alpha = deg2rad(config.alpha_deg);
    Vinf_vec = config.V_inf * [cos(alpha) 0 sin(alpha)];
    Vinf_hat = Vinf_vec / norm(Vinf_vec);
    Lref = 1e4 * b;

    boundA = collect_field(panels, 'boundA');
    boundB = collect_field(panels, 'boundB');
    A_inf  = boundA + Lref*Vinf_hat;
    B_inf  = boundB + Lref*Vinf_hat;

    xg = linspace(-1.5*croot, 2.5*croot, 6);
    yg = linspace(-0.95*b/2, 0.95*b/2, 11);
    zg = linspace(-0.25*b, 0.25*b, 5);
    [X, Y, Z] = meshgrid(xg, yg, zg);

    U = Vinf_vec(1) * ones(size(X));
    V = Vinf_vec(2) * ones(size(X));
    W = Vinf_vec(3) * ones(size(X));

    pts = [X(:), Y(:), Z(:)];
    for k = 1:size(pts,1)
        Vtot = biot_savart_vec(pts(k,:), A_inf, boundA) + ...
               biot_savart_vec(pts(k,:), boundA, boundB) + ...
               biot_savart_vec(pts(k,:), boundB, B_inf);
        Vind = sum(Gamma .* Vtot, 1);
        U(k) = U(k) + Vind(1);
        V(k) = V(k) + Vind(2);
        W(k) = W(k) + Vind(3);
    end
end

% ------------------------------------------------------------------------
%  VISUALIZATION DASHBOARD
% ------------------------------------------------------------------------
function plot_dashboard(panels, Gamma, results, config)
    N = numel(panels);
    fig = figure('Position', [100, 100, 1400, 700], 'Color', 'w');

    % Build a single Faces/Vertices patch representation shared by both
    % 3D subplots.
    Vtx = zeros(4*N, 3);
    Fcs = zeros(N, 4);
    for i = 1:N
        i0 = (i-1)*4;
        Vtx(i0+1,:) = panels(i).P1;
        Vtx(i0+2,:) = panels(i).P2;
        Vtx(i0+3,:) = panels(i).P3;
        Vtx(i0+4,:) = panels(i).P4;
        Fcs(i,:) = i0 + (1:4);
    end

    %% ---- Subplot 1: 3D surface pressure distribution ----
    ax1 = subplot(1,3,1);
    patch(ax1, 'Faces', Fcs, 'Vertices', Vtx, 'FaceVertexCData', results.dCp, ...
          'FaceColor', 'flat', 'EdgeColor', [0.25 0.25 0.25], 'EdgeAlpha', 0.35);
    colormap(ax1, turbo);
    cb = colorbar(ax1);
    cb.Label.String = '\Delta C_p';
    axis(ax1, 'equal'); grid(ax1, 'on'); view(ax1, -35, 25);
    xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)'); zlabel(ax1,'z (m)');
    title(ax1, sprintf('%s Wing - \\Delta C_p Distribution', config.wing_type), ...
          'Interpreter', 'tex');
    camlight(ax1, 'headlight'); lighting(ax1, 'gouraud'); material(ax1, 'dull');

    %% ---- Subplot 2: 3D streamlines / wingtip vortex roll-up ----
    ax2 = subplot(1,3,2);
    hold(ax2, 'on');
    patch(ax2, 'Faces', Fcs, 'Vertices', Vtx, 'FaceColor', [0.78 0.78 0.85], ...
          'EdgeColor', [0.35 0.35 0.35], 'FaceAlpha', 0.65);

    [X, Y, Z, U, V, W] = build_flow_grid(panels, Gamma, config);
    [sy, sz] = meshgrid(linspace(-0.9*config.span/2, 0.9*config.span/2, 9), ...
                         linspace(-0.05*config.span, 0.05*config.span, 3));
    sx = (-1.4*config.root_chord) * ones(size(sy));
    hlines = streamline(ax2, X, Y, Z, U, V, W, sx(:), sy(:), sz(:));
    set(hlines, 'Color', [0.10 0.40 0.90], 'LineWidth', 1.0);

    axis(ax2, 'equal'); grid(ax2, 'on'); view(ax2, -35, 25);
    xlabel(ax2,'x (m)'); ylabel(ax2,'y (m)'); zlabel(ax2,'z (m)');
    title(ax2, 'Streamlines & Wingtip Vortex Roll-up');
    hold(ax2, 'off');

    %% ---- Subplot 3: Performance readout + spanwise lift distribution ----
    ax3 = subplot(1,3,3);
    axis(ax3, 'off');
    txt = { ...
        sprintf('PLANFORM: %s', config.wing_type), ...
        sprintf('Span b = %.2f m', config.span), ...
        sprintf('Root Chord = %.2f m', config.root_chord), ...
        sprintf('Tip Chord = %.2f m', config.tip_chord), ...
        sprintf('S_{ref} = %.3f m^2', results.S_ref), ...
        '', ...
        sprintf('\\alpha = %.1f deg,   V_\\infty = %.1f m/s', config.alpha_deg, config.V_inf), ...
        sprintf('Dynamic Pressure q = %.1f Pa', results.q_inf), ...
        '', ...
        sprintf('Lift L = %.1f N', results.L_total), ...
        sprintf('C_L = %.4f', results.CL), ...
        sprintf('C_{Di} = %.5f', results.CDi), ...
        sprintf('L/D = %.2f', results.LD) };
    text(ax3, 0.02, 0.98, txt, 'VerticalAlignment', 'top', 'FontSize', 11, ...
         'FontName', 'Consolas', 'Interpreter', 'tex');

    axInset = axes('Position', [0.70 0.12 0.27 0.32]);
    plot(axInset, results.y_span, results.lift_per_unit_span, '-o', ...
         'LineWidth', 1.5, 'MarkerSize', 3, 'Color', [0.85 0.20 0.20]);
    grid(axInset, 'on');
    xlabel(axInset, 'Span y (m)'); ylabel(axInset, 'Lift / span (N/m)');
    title(axInset, 'Spanwise Lift Distribution', 'FontSize', 9);

    sgtitle(fig, sprintf('3D Vortex Lattice Method - %s Wing Aerodynamic Analysis', ...
            config.wing_type));
end