%s = DESIGN(varargin)
%       Transfer model of extractor to design its geometry and other operating parameters for given influent and effluent compositions
%       This function first calls one of the subroutines (depending on extractor type) to calculate parameters of mass transfer
%       Substance properties dependent on temperature are calculated using getProperties()
%       First, the calculation of transfer model without axial dispersion is performed (assuming plug flow)
%       Using the result of previous calculation as an initial estimate, a model with axial dispersion is calculated
%       Calculation method is described in "Design of Liquid-Liquid Extractors" (1997) by S. Hartland and A. Kumar
%
%PARAMETERS:
%       There are many optional parameters that can be modified, otherwise calculation will proceed with their default values
%       See documentation for list of all possible parameters as well as their default values
%       Parameters can be entered as field & value pairs:
%              design('param1Name', param1Value, 'param2Name', param2Value, ...)
%       Alternatively, you can input a pre-prepared struct as the FIRST argument (to save time with entering constant parameters all over again)
%              design(STRUCT, 'param1Name', param1Value, ...)
%
%OUTPUT:
%       Is the struct object with initial data as well as calculation results
%       See documentation for list of all result properties
%
%CALL EXAMPLES:
% s=design('verb',2, 'type','packed', 'phase1',{'SA40',26500,0.93e-2,1e-4}, 'phase2',{'TOL',6000,0,NaN}, 't',60, 'packing', struct('type','Raschig', 'Dn',51e-3, 'at',100, 'eps',0.94));
% s=design('verb',1, 'type','packed', 'L',3, 'phase1',{'SA40',26500,0.93e-2,NaN}, 'phase2',{'TOL',6000,0,NaN}, 't',60, 'packing',12);
% s=design('verb',1, 'type','rotDisc', 'dispPhase',2, 'phase1',{'SA40',26500,0.93e-2,1e-4}, 'phase2',{'TOL',6000,0,NaN}, 't',60, 'N',0.2);
% s=design('verb',1, 'type','karr', 'phase1',{'SA40',26500,0.93e-2,1e-4}, 'phase2',{'TOL',6000,0,NaN}, 't',60, 'N',2.44, 'stroke',8.20e-3, 'H',50e-3);


function s = design(varargin)
    timestamp = now;
    addpath('lib'); % dependencies

    s = extendStruct(varargin{:}); % initiate 's', which is a struct with all input properties, serving as the universal state object
    s = designPrepareState(s); % see designPrepareState() definition below

    % OVERRIDES
    % interfacial tension
    doTensionCalc = isfield(s, 'tensionF');
    ~doTensionCalc && s.verb >= 2 && fprintf(1, 'OVERRIDDEN gamma = %.3f N/m (interfacial tension)\n', s.tension);
    % column diameter
    if isfield(s, 'Dt'); DtOverride = s.Dt; else; DtOverride = 0; end
    DtOverride > 0 && s.verb >= 2 && fprintf(1, 'OVERRIDDEN Dt = %.3f m (column diameter)\n', s.Dt);
    % column length OVERRIDE
    if isfield(s, 'L'); LOverride = s.L; else; LOverride = 0; end
    LOverride > 0 && s.verb >= 2 && fprintf(1, 'OVERRIDDEN L = %.3f m (column height)\n', s.L);

    % whether mass flows are unknown
    mFlowUnknown = isnan(s.phase1.m) || isnan(s.phase2.m);

    % keys of 's' that point to s.phase1 and s.phase2
    dispKey = ['phase' num2str(  s.dispPhase)];
    contKey = ['phase' num2str(3-s.dispPhase)];

    % obtain appropriate equilibrium function as u(w), where 'w' is CPL mass frac in aq, 'u' in org
    sol = s.phase2.id; % solvent as string
    SA  = str2num(s.phase1.id(3:4)); % SA percentage as number
    w_aq = [1e-4:1e-4:0.01, 0.015:5e-3:0.1]; w_aq = [-w_aq(end:-1:1) 0 w_aq]; % tabelate in a very wide range, but with decreasing density
    [w_aq, w_org, w_aq_max, w_org_max] = getEqFun(sol, SA, s.t, w_aq); % get eqFun as tabelated function

    % extrapolate the last points to infinity and beyond!
    w_aq = [-1e4; w_aq; 1e4];
    w_org = [w_org(1)*1e5; w_org; w_org(end)*1e5];
    % the points needs to exist, because bvp solver sometimes looks there and needs at least something, so it can safely return to the sane region

    s.verb >= 2 && headline('PHYSICAL PROPERTIES');
    msg = ['Continuous phase: %s\nDispersed phase: %s\nt = %.1f °C\n'...
        'rhoC = %.0f kg/m3 (density)\nrhoD = %.0f kg/m3\nmuC = %.2e Pa*s (dynamic viscosity)\nmuD = %.2e Pa*s\nDC = %.2e m2/s (CPL diffusivity)\nDD = %.2e m2/s\n\n'];
    s.verb >= 2 && fprintf(1, msg, s.(contKey).id, s.(dispKey).id, s.t, s.rhoC, s.rhoD, s.muC, s.muD, s.DC, s.DD);

    s.verb >= 2 && headline('PRE-RESULTS');

    % mass balance is calculated now, unless length is overriden, in which case it is calculated after the transfer model
    w_warning_issued = false;
    function callStage(str)
        % call the stage model to complete mass balance, get minimal solvent flow and the exact NTP value
        stageResults = extractor([s.phase1.m, s.phase2.m, s.phase1.wi, s.phase1.we, s.phase2.wi, s.phase2.we], sol, SA, s.t);
        s.phase1.m = stageResults.mF; s.phase2.m = stageResults.mS;
        wR = stageResults.wR; wF = stageResults.wF;
        wS = stageResults.uS; wE = stageResults.uE;

        % warning if concentrations are out of data bounds
        if ~w_warning_issued
            msg = 'EXTRAPOLATION\nrequested CPL mass fraction in %s phase (%.2f %%) exceeds the available data (%.2f %%)';
            warning('off', 'backtrace');
            if wF > w_aq_max  && s.verb >= 1; warning(msg, 'aqueous', 100*wF, 100*w_aq_max ); end
            if wE > w_org_max && s.verb >= 1; warning(msg, 'organic', 100*wE, 100*w_org_max); end
            warning('on', 'backtrace');
            w_warning_issued = true;
        end

        msg = '%smass balance:\n   wF = %.4f %% (feed)\n   wR = %.4f %% (raffinate)\n   wS = %.4f %% (solvent)\n   wE = %.4f %% (extract)\n';
        s.verb >= 1 && fprintf(1, msg, str, wF*100, wR*100, wS*100, wE*100);
        s.verb >= 1 && mFlowUnknown && fprintf(1, '   m1 = %.2e kg/h (aqueous phase mass flow)\n   m2 = %.2e kg/h (organic phase mass flow)\n', s.phase1.m, s.phase2.m);
        s.verb >= 1 && fprintf(1, 'mSmin = %.3e kg/hr\n', stageResults.mSmin);
        s.verb >= 2 && fprintf(1, '   actual flow is %.0f %% of mSmin\n', s.phase2.m / stageResults.mSmin * 100);
        s.verb >= 1 && fprintf(1, 'NTP = %.2f (number of theoretical plates)\n\n', stageResults.TP);
    end

    if LOverride == 0; callStage(''); % call mass balance
    else wR = s.phase1.we; wF = s.phase1.wi; wS = s.phase2.wi; wE = s.phase2.we; end % or just fill in fallback values
    if isnan(wR); wR = 0; end
    if isnan(wF); wF = 1e-2; end
    if isnan(wS); wS = 0; end
    if isnan(wE); wE = 1e-2; end

    s.QD = s.(dispKey).m / 3600 / s.rhoD; % volume flow of disp [m3/s]
    s.QC = s.(contKey).m / 3600 / s.rhoC; % volume flow of cont [m3/s]
    % note: phase mass flows are not converted from [kg/hr], because they won't be needed afterwards
    s.R = s.QD / s.QC; % disp / cont ratio of superficial velocities


    %{
        CHOOSE CALCULATION SUBPROGRAM FOR THE SELECTED INSTRUMENT
        Following s fields will be calculated:
            aint (interfacial surface [m2/m3]), A (column cross-sectional area [m2]),
            DAD, DAC (dispersion coeficients [m2/s]), betaD, betaC (mass transfer coeficients [m/s])
            phi (dispersed phase retention, only for the "check for dispersion" part)
        Now we will assign insfun(), a "pointer" to the appropriate instrument function
    %}
    fields = {'packed', 'rotDisc', 'kuhni', 'karr'}; % , 'packed2'
    if     strcmpi(s.type, fields{1}); insfun = @(s,dt) designPacked (s,dt);
    elseif strcmpi(s.type, fields{2}); insfun = @(s,dt) designRotDisc(s,dt);
    elseif strcmpi(s.type, fields{3}); insfun = @(s,dt) designRotDisc(s,dt);
    elseif strcmpi(s.type, fields{4}); insfun = @(s,dt) designKarr   (s,dt);
    % elseif strcmpi(s.type, fields{5}); insfun = @(s,dt) designPacked2 (s,dt);
    else error('column type ''%s'' is not implemented. Permitted values of ''type'' are (case insensitive):\n%s', s.type, strjoin(fields, ', '));
    end


    %{
        CALL INSFUN
        first in sizing mode, to determine Dt and VCE
        then in design mode to calculate aint, DAD, DAC, betaD, betaC and phi
    %}
    % get mean tension, maximal Dt and minimal VCE (all of them throughout whole column)
    if doTensionCalc
        w = (wR : 1e-4 : wF); % tabelate aqueous mass frac
        tensionVec = zeros(size(w)); % vector of interfacial tensions in each section [N/m]
        maxDt = 0; % maximal column diameter [m]
        minVCE = Inf; % minimal flooding velocity (cont) [m/s]
        for i = 1:length(w)
            s.tension = s.tensionF(w(i)); % calculate current tension based on aq mass frac
            tensionVec(i) = s.tension;
            s = insfun(s, true); % call instrument function to calculate Dt and VCE of current section
            if s.Dt  > maxDt;  maxDt  = s.Dt;  end % save current value of Dt so as to find the maximum
            if s.VCE < minVCE; minVCE = s.VCE; end % minimum VCE (analogous)
        end
        s.Dt = maxDt;
        s.VCE = minVCE;
        s.tension = mean(tensionVec);
        s.verb >= 2 && fprintf(1, 'gamma = %.3f N/m (CALCULATED interfacial tension)\n', s.tension);
    % tension is already set, so call instrument function to calculate Dt and VCE for the whole column
    else
        s = insfun(s, true);
    end

    if DtOverride > 0; s.Dt = DtOverride; % value of Dt that was set as a parameter is returned back (the previous calculations were necessary to obtain VCE as informative value)
    else s.verb >= 1 && fprintf(1, 'Dt = %.3f m (CALCULATED column diameter)\n', s.Dt); end

    A = pi * s.Dt^2 / 4;
    VC = s.QC / A;
    VD = s.R * VC;
    msg = 'QCE = %.2e m3/s (cont phase FLOODING flow at safety factor %.0f%%)\nQC = %.2e m3/s (volume flow)\nQD = %.2e m3/s\nR = %.2f (flow ratio = QD/QC)\n\n';
    s.verb >= 1 && fprintf(1, msg, s.VCE*A, VC/s.VCE*100, s.QC, s.QD, s.R);
    msg = 'VCE = %.2e m/s (cont phase FLOODING velocity)\nVC =  %.2e m/s (velocity of cont phase)\nVD =  %.2e m/s (velocity of disp phase)\n\n';
    s.verb >= 2 && fprintf(1, msg, s.VCE, VC, VD);
    if VC > s.VCE; warning('off', 'backtrace'); warning('Column is predicted to be in FLOODING!'); warning('on', 'backtrace'); end

    s = insfun(s, false); % call instrument function only once to calculate: aint, A, DAD, DAC, betaD, betaC, phi
    Ksi = s.phi/(1 - s.phi) * (s.rhoD / s.rhoC * s.muC / s.muD)^0.3; % dimensionless number to check for dispersion stability (according to Perry), evaluation is in RESULTS section


    %{
        THE MAIN PART OF THE PROGRAM - extractor transfer model
        Calculated using Matlab BVP solver on the differential equations
        First, all the necessary nested functions are defined
        Then, the plug flow model is calculated (code labeled PLUGFLOW)
        Finally, the model with axial dispersion is calculated using plugflow results as estimate (code labeled DISPERSION)
        Key for variable names:     x ~ aq ~ phase1     y ~ org ~ phase2
    %}
    % assign calculated mass transfer variables
    aux = {'D', 'C'};
    aqKey = aux{s.dispPhase}; orgKey = aux{3-s.dispPhase};
    Qaq  = s.(['Q' aqKey]);
    Qorg = s.(['Q' orgKey]);
    DAaq  = s.(['DA' aqKey]);
    DAorg = s.(['DA' orgKey]);
    rhoaq  = s.(['rho' aqKey]);
    rhoorg = s.(['rho' orgKey]);
    betaaq  = s.(['beta' aqKey]);
    betaorg = s.(['beta' orgKey]);

    % assign concentration boundaries
    cx1 = wR * rhoaq;  cx2 = wF * rhoaq;
    cy1 = wS * rhoorg; cy2 = wE * rhoorg;

    % get data pairs [cx1, cx2], [cy1, cy2] around specified 'cx' OR 'cy' value, the other must be specified as NaN
    function [cxx, cyy] = get_data_pairs(cx, cy)
        if isnan(cy)
            wx_curr = cx / rhoaq;
            if wx_curr > w_aq(end) || wx_curr < w_aq(1); error('Internal error\nwx_curr = %.3e\n', wx_curr); end
            j = find(w_aq > wx_curr);
        else
            wy_curr = cy / rhoorg;
            if wy_curr > w_org(end) || wy_curr < w_org(1); error('Internal error\nwy_curr = %.3e\n', wy_curr); end
            j = find(w_org > wy_curr);
        end
        j = j(1) - 1; % index (in w_aq and w_org) that points to the lower w_aq-w_org pair
        if j < 1; j = 1; end
        cxx = rhoaq  * [w_aq(j),  w_aq(j+1)];
        cyy = rhoorg * [w_org(j), w_org(j+1)];
    end

    % get m & q (slope & intercept) of tangent line to the eq curve cy(cx) at current cx, cy
    function [m, q] = get_mq(cx)
        [cxx, cyy] = get_data_pairs(cx, NaN);
        m = diff(cyy) / diff(cxx);
        q = cyy(1) - m * cxx(1);
    end

    % equilibrium function   cy(cx)
    function cy = eqFun(cx)
        [cxx, cyy] = get_data_pairs(cx, NaN);
        cy = cyy(1) + diff(cyy) * (cx - cxx(1)) / diff(cxx);
    end
    % equilibrium function inverse   cx(cy)
    function cx = eqFunInv(cy)
        [cxx, cyy] = get_data_pairs(NaN, cy);
        cx = cxx(1) + diff(cxx) * (cy - cyy(1)) / diff(cyy);
    end

    % callback for PLUGFLOW (1st order differential equation)
    function dc = diffEq1(z, c, L)
        cx = c(1); cy = c(2);

        [m, q] = get_mq(cx);
        cyeq = eqFunInv(cy);
        Kx = 1/((1/betaaq) + (1/betaorg/m)); % coefficient of mass transfer [m/s]
        Tx = L * Kx * s.aint * s.A / Qaq;

        dc = zeros(2, 1);
        dc(1) = -Tx * (cx - cyeq);
        dc(2) = -Tx * (cx - cyeq) * Qaq / Qorg;
    end

    % callback for DISPERSION (2nd order differential equation)
    function dc = diffEq2(z, c, L)
        cx = c(1); cy = c(2); dcx = c(3); dcy = c(4);

        [m, q] = get_mq(cx);
        cyeq = eqFunInv(cy);
        Kx = 1/((1/betaaq) + (1/betaorg/m)); % coefficient of mass transfer [m/s]
        Tx = L * Kx * s.aint * s.A / Qaq;

        dc = zeros(4, 1);
        dc(1) = dcx;
        dc(2) = dcy;
        dc(3) =  Pex(L)*dcx + Pex(L)*Tx*(cx - cyeq);
        dc(4) = -Pey(L)*dcy - Pey(L)*Tx*(cx - cyeq) * Qaq / Qorg;
    end
    Pex = @(L) L * Qaq  / s.A / DAaq;
    Pey = @(L) L * Qorg / s.A / DAorg;


    %{
        PLUGFLOW
        z = dimensionless length
        c = [cx, cy]
    %}
    % callback with guess solution as linear 'c' profiles throughout 'z'
    guess = @(z) [
        cx2 - z*(cx2-cx1)
        cy2 - z*(cy2-cy1)];

    % call BVP solver
    % IF calculate length, ELSE calculate raffinate with defined length
    if LOverride == 0
        % callback with boundary conditions
        boundconds = @(ca,cb,L) [
            ca(1) - cx2
            cb(1) - cx1
            cb(2) - cy1];

        L0 = stageResults.TP * 0.5; % initial L using a very crude estimate of HETP = 0.5 m
        solinit = bvpinit(0:0.01:1, guess, L0);
        bvpsol1 = bvp5c(@diffEq1, boundconds, solinit);
        s.L = bvpsol1.parameters;
    else
        % no raffinate definition when length is overriden
        boundconds = @(ca,cb) [
            ca(1) - cx2
            cb(2) - cy1];

        solinit = bvpinit(0:0.01:1, guess);
        bvpsol1 = bvp5c(@(z,c) diffEq1(z, c, s.L), boundconds, solinit);
    end


    %{
        DISPERSION
        z = dimensionless length
        c = [cx, cy, dcx, dcy]
    %}
    options = bvpset('NMax', 10000);

    % callback with guess solution using the plugflow solution
    function c = guess2(z)
        % find index that points to the lower number
        j = find(bvpsol1.x > z);
        if isempty(j) || j(1) == length(bvpsol1.x); j = length(bvpsol1.x) - 1;
        elseif j(1) == 1; j = 2;
        else j = j(1) - 1; end
        c1 = bvpsol1.y(:, j); c2 = bvpsol1.y(:, j+1);
        z1 = bvpsol1.x(j); z2 = bvpsol1.x(j+1);

        % interpolation within the original points
        c = c1 + (c2-c1) .* (z - z1) / (z2-z1);

        % and a simple two point derivation
        dc = (bvpsol1.y(:, j+1) - bvpsol1.y(:, j)) / (z2-z1);
        c = [c; dc];
        % old guess - linear profiles and their derivatives
        % c = cx2 - z*(cx2-cx1); cy2 - z*(cy2-cy1); -(cx2-cx1); -(cy2-cy1)];
    end

    % call BVP solver
    % IF calculate length, ELSE calculate raffinate with defined length
    if LOverride == 0
        % callback with boundary conditions
        boundconds2 = @(ca,cb,L) [
            cb(1) - cx1
            ca(3) - Pex(L) * (ca(1) - cx2)
            cb(3)
            ca(4)
            cb(4) - Pey(L) * (cy1 - cb(2))];

        solinit = bvpinit(0:0.01:1, @guess2, s.L * s.guessFactor); % see default.guessFactor
        bvpsol2 = bvp5c(@diffEq2, boundconds2, solinit, options);
        s.L2 = bvpsol2.parameters;
    else
        % different boundary conds for column length override
        boundconds2 = @(ca,cb) [
            ca(3) - Pex(s.L) * (ca(1) - cx2)
            cb(3)
            ca(4)
            cb(4) - Pey(s.L) * (cy1 - cb(2))];

        solinit = bvpinit(0:0.01:1, @guess2);
        bvpsol2 = bvp5c(@(z,c) diffEq2(z, c, s.L), boundconds2, solinit, options);
    end
% DEBUG
% fprintf(1, 'PLUGFLOW RESIDUAL   = %.3e\n', abs(Qaq * (bvpsol1.y(1,1) - bvpsol1.y(1,end)) - Qorg * (bvpsol1.y(2,1) - bvpsol1.y(2,end))));
% fprintf(1, 'DISPERSION RESIDUAL = %.3e\n', abs(Qaq * (bvpsol2.y(1,1) - bvpsol2.y(1,end)) - Qorg * (bvpsol2.y(2,1) - bvpsol2.y(2,end))));
% yy2=bvpsol2.y*0; for i = 1:length(bvpsol2.x); yy2(:, i) = guess2(bvpsol2.x(i)); end
% figure();plot(bvpsol2.x, bvpsol2.y, '-', bvpsol2.x, yy2,'--', [0,0,1,1],[cx2,cy2,cx1,cy1],'xk');grid on;
% {'','MB','plug','disp';'R',cx1,bvpsol1.y(1,end),bvpsol2.y(1,end);'F',cx2,bvpsol1.y(1,1),bvpsol2.y(1,1);'S',cy1,bvpsol1.y(2,end),bvpsol2.y(2,end);'E',cy2,bvpsol1.y(2,1),bvpsol2.y(2,1)}


    %{
      FINISH THE PROGRAM
    %}
    s.verb >= 2 && headline('FINAL RESULTS');

    % evaluate dispersion stability (according to Perry), calculation is in the CALL INSFUN section
    if     Ksi < 0.3; msg = 'confirmed (Ksi < 0.3)';
    elseif Ksi < 0.5; msg = 'probably confirmed (Ksi 0.3 - 0.5)';
    elseif Ksi < 2.0; msg = 'neither confirmed nor disproven (Ksi 0.5 - 2.0)';
    elseif Ksi < 3.3; msg = 'probably disproven (Ksi 2.0 - 3.3)';
    else              msg = 'disproven (Ksi > 3.3)'; end

    s.verb >= 2 && fprintf(1, 'maximal Ksi = %.2f, therefore definition of %s as dispersed phase is %s\n\n', Ksi, s.(dispKey).id, msg);

    % calculate number of mixers and total power
    if isfield(s, 'H') && isfield(s, 'P')
        L = s.L; if isfield(s, 'L2') && ~isnan(s.L2) && s.L2 > L; L = s.L2; end
        s.mixers = L / s.H;
        s.Ptotal = ceil(s.mixers) * s.P;
        s.verb >= 1 && fprintf(1, 'number of mixers = %.1f, rounded to %d\ntotal power P = %.0f W\n\n', s.mixers, ceil(s.mixers), s.Ptotal);
    end

    if LOverride == 0;
        s.verb >= 1 && fprintf(1, 'L = %.3f m (column height - PLUGFLOW model)\n', s.L);
        s.verb >= 2 && fprintf(1, '    HETP = %.3f (height equivalent of theoretical plate)\n', s.L / stageResults.TP);
        s.verb >= 1 && fprintf(1, 'L = %.3f m (column height - DISPERSION model)\n', s.L2);
        s.verb >= 2 && fprintf(1, '    HETP = %.3f (height equivalent of theoretical plate)\n', s.L2 / stageResults.TP);
        s.verb >= 1 && fprintf(1, '\n');
    else
        % fill in the calculate raffinate and complete mass balance
        s.phase1.we = bvpsol1.y(1,end) / rhoaq;
        callStage('PLUGFLOW MODEL ');
        s.phase1.we = bvpsol2.y(1,end) / rhoaq;
        callStage('DISPERSION MODEL ');
    end
    s.verb >= 2 && fprintf(1, 'Execution of design.m took %.2f s\n\n', (now-timestamp)*24*3600);
end


% validates & completes input struct with default values and returns it
function s = designPrepareState(s)
    % required fields
    if ~isfield(s, 'type'); error('Parameter ''type'' is required (string definition of column type)'); end
    if ~isfield(s, 'phase1'); error('Parameter ''phase1'' is required (as a struct which defines the feed phase)'); end
    if ~isfield(s, 'phase2'); error('Parameter ''phase2'' is required (as a struct which defines the solvent phase)'); end

    % make phase struct from phase cell
    function ps = phaseStruct(pc)
        if numel(pc) ~= 4; error('When phase is defined by a cell array, it must have 4 values: {phase id, mass flow [kg/hr], mass frac in, mass frac out}'); end
        ps = struct('id',pc{1}, 'm',pc{2}, 'wi',pc{3}, 'we',pc{4});
    end
    if iscell(s.phase1); s.phase1 = phaseStruct(s.phase1); end
    if iscell(s.phase2); s.phase2 = phaseStruct(s.phase2); end

    % struct field of default input values
    default = struct();
    default.verb = 1; % how much console output will there be
    default.flowSafety = 0.50; % safety factor for flooding
    default.dispPhase = 1; % which phase is dispersed (either 1 or 2)
    default.t = 60; % temperature [°C]

    % multiplier of guess column height for axial dispersion model in design mode
    % change it when solution is unstable (no sensible solution is obtained) or when computation time is too long
    default.guessFactor = 2;

    % defaults for packed columns - define packing properties
    if strcmpi(s.type, 'packed') || strcmpi(s.type, 'packed2')
        % default value - predefined packing
        if ~isfield(s, 'packing'); s.packing = 5; end
        % full definition of packing received, check whether all fields are in place:
        if isstruct(s.packing)
            fields = {'type', 'Dn', 'at', 'eps'};
            for i = 1:numel(fields)
                if ~isfield(s.packing, fields{i}); error('Packing must have all of the following properties defined: %s', strjoin(fields, ', ')); end
            end
        % or get values from a predefined table of packings using packing.m:
        else
            res = packing(s.packing);
            s.packing = struct('type',res{2}, 'Dn',res{3}/1000, 'at',res{6}, 'eps',res{7}/100);
        end
    % defaults for rotating disc extractor
    elseif strcmpi(s.type, 'rotDisc')
        default.N = 0.3; % frequency of rotation [Hz]
    % defaults for Karr extractor
    elseif strcmpi(s.type, 'Karr')
        default.N = 4; % frequency of strokes [Hz]
        default.stroke = 15e-3; % length of strokes [m]
        default.aperture = 0.58; % fractional free area of discs
        default.H = 50e-3; % compartment height [m]
    end

    % complete state object with default values if any of the fields are missing
    fields = fieldnames(default);
    for i = 1:numel(fields)
        if ~isfield(s, fields{i})
            s.(fields{i}) = default.(fields{i});
        end
    end

    % validate & complete structs with phase definitions
    if s.dispPhase ~= 1 && s.dispPhase ~= 2; error('''dispPhase'' must be either 1 or 2'); end
    if ~isfield(s.phase1, 'id') || ~isfield(s.phase2, 'id'); error('Both phases must have ''id'' property defined'); end
    if strcmpi(s.phase1.id, s.phase2.id); error('Phases must be different, both are ''%s''', s.phase1.id); end
    fields = {'m', 'wi', 'we'}; % mass balance quantities
    m = 0; % number of missing values - exactly one mass balance quantity should be missing or NaN
    for i = 1:numel(fields)
        if ~isfield(s.phase1, fields{i}); s.phase1.(fields{i}) = NaN; end
        if ~isfield(s.phase2, fields{i}); s.phase2.(fields{i}) = NaN; end
        if isnan(s.phase1.(fields{i})); m = m + 1; end
        if isnan(s.phase2.(fields{i})); m = m + 1; end
    end
    n = isfield(s, 'L')+1; % one MB quantity is unknown if length is defined, two if length is defined
    if m ~= n; error('Both phases must have following fields defined: %s\nof which exactly %d must be NaN, but %d NaN values encountered', strjoin(fields, ', '), n, m); end
    if isfield(s, 'L') && (isnan(s.phase1.m) || isnan(s.phase2.m)); error('When length is overriden, all mass flows must be specified (while 2 of the mass frac values must be NaN)'); end


    % assign physical properties using correlations on data in 'data/properties.csv'
    IDD = s.(['phase' num2str(  s.dispPhase)]).id; % id disp phase
    IDC = s.(['phase' num2str(3-s.dispPhase)]).id; % id cont phase
    data = readCSVasCell('data/properties.csv');
    getRow = @(label) [data{find(strcmpi(data(:, 1), label)), 3:end}]; % find data row by its label, and return the rest of the row as vector

    % assign temperature dependent phase properties: density [kg/m3], viscosity [Pa*s], diffusivity [m2/s]
    % they are considered to be constant throughout the whole calculation - since they are bulk properties, they should be independent on w(CPL)
    fields = {'rho', 'mu', 'D'};
    tvals = getRow('t');
    for i = 1:numel(fields)
        Dvals = getRow([fields{i} '_' IDD]);
        Cvals = getRow([fields{i} '_' IDC]);
        Dfit = polyfit(tvals, Dvals, 2);
        Cfit = polyfit(tvals, Cvals, 2);
        s.([fields{i} 'D']) = polyval(Dfit, s.t);
        s.([fields{i} 'C']) = polyval(Cfit, s.t);
    end

    if isfield(s, 'gamma'); s.tension = s.gamma; end % clear up ambiguity of paramater name
    if isfield(s, 'tension'); return; end
    % create function for interfacial tension
    label = [s.phase1.id '_' s.phase2.id];
    label = regexprep(label, 'SA\d{2}', 'SA40'); % a kind of hotfix for missing SA42 data
    tensionFit = polyfit(getRow(['w_CPL_' label]), getRow(['gamma_' label]), 2);
    s.tensionF = @(w) polyval(tensionFit, w);
end
