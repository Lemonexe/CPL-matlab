%s = EXTRACTOR(v, sol, SA, t)
%       calculates mass balance of extractor and simulation of theoretical plates using the McCabe-Thiele method
%       output is a struct object containing all kinds of data, most importantly the number of theoretical plates (TP)
%
%PARAMETERS:
%       v = [mF, mS, wF, wR, uS, uE]
%              with 'm' being mass or mass flow, 'w' and 'u' are mass fractions in aqueous and organic phases
%              F = feed, S = solvent, R = raffinate, E = extract
%              it is absolutely imperative that any member of vector 'v', but exactly one, is NaN
%              that is the variable which will be calculated from mass balance
%       sol = organic solvent = 'TOL', 'TCE'
%       SA = content of ammonium sulphate % = 30, 35, 40
%       t = temperature in °C, should be between 30 and 60, but extrapolation beyond these bounds is possible
%
%CALL EXAMPLE:
%       extractor([100, 10, 0.01, 0.001, 0, NaN], 'TOL', 40, 60)

function s = extractor(v, sol, SA, t)
    % s = structured array with everything - state data & results
    s = struct('v', v, 'sol', sol, 'SA', SA, 't', t);

    % utility functions - convert relative and absolute fractions
    s.abs2rel = @(w) w./(1 - w);
    s.rel2abs = @(W) W./(1 + W);

    %{
    1. calculate mass balance
    2. check if the (now completed) balance data makes sense
    3. get tabelated equilibrium function for given temperature (see getEqFun.m)
    4. calculate minimal solvent flow
    5. McCabe-Thiele
        see the function definitions below for more details
    %}
    s = balance(s);
    checkConstraints(s);
    [s.aq, s.org] = getEqFun(s.sol, s.SA, s.t);
    s = minFlow(s);
    s = McCabeThiele(s);
end

% perform initial mass balance
function s = balance(s)
    % at this point we need to determine which one of the v members is omitted
    s.i = find(isnan(s.v));
    if isempty(s.i) % too many inputs
        error('Balance is overdetermined - one input value must be omitted!');
    elseif length(s.i) > 1 % not enough inputs
        error('Balance is incomplete - there must be 5 input values!');
    end

    % fsolve for mR, mE and the one missing input
    s.v(s.i) = fzero(@balanceCallback, 0.5, [], s.v, s.i);

    % unpack v (now completed with fsolved value) into named variables
    v = num2cell(s.v);
    [s.mF, s.mS, s.wF, s.wR, s.uS, s.uE] = v{:};

    s.mR = s.mF * (1-s.wF) / (1-s.wR); % aq balance
    s.mE = s.mS * (1-s.uS) / (1-s.uE); % org balance

    % from absolute fractions to relative
    s.ma = s.mF*(1 - s.wF);
    s.mt = s.mS*(1 - s.uS);
    s.WF = s.abs2rel(s.wF);
    s.WR = s.abs2rel(s.wR);
    s.US = s.abs2rel(s.uS);
    s.UE = s.abs2rel(s.uE);
end

% fsolve callback (residual function). Also needs 'v' vector and 'i', the index of missing input
function res = balanceCallback(x, v, i)
    % unpack v (now completed with temporary value) into named variables
    v(i) = x;
    v = num2cell(v);
    [mF, mS, wF, wR, uS, uE] = v{:};

    mR = mF * (1-wF) / (1-wR); % aq balance
    mE = mS * (1-uS) / (1-uE); % org balance
    res = mF*wF + mS*uS - mR*wR - mE*uE; % residue for CPL balance
end

% check consistency of the mass balance
function checkConstraints(s)
    if     s.wR > s.wF; error('wR > wF \n%.4g > %.4g\n(raffinate has to have lower concentration than feed)',  s.wR, s.wF);
    elseif s.uS > s.uE; error('uS > uE \n%.4g > %.4g\n(solvent has to have lower concentration than extract)', s.uS, s.uE);
    end

    % check if any of balance values < 0
    varNames = {'mF', 'mS', 'wF', 'wR', 'uS', 'uE'};
    for i = 1:length(varNames)
        name = varNames{i};
        val = s.(name);
        if val < 0; error('NEGATIVE value encountered: %s = %.4g', name, val); end
    end
end

% implementation of bisection method to find minimal possible 'mS', which means finding lowest 'mt' where operation line intersects equilibrium line
function s = minFlow(s)
    mt1 = s.ma*10; mt2 = 0; % upper and lower limit of 'mt' in current iteration

    if s.wR > max(s.aq)
        warning('Minimal mass of solvent cannot be calculated, as wR exceeds experimental data range for aqueous phase.\nwR = %.4g\nmax w(aq) = %.4g', s.wR, max(s.aq));
        s.mSmin = 0; return
    end

    % difference between operation line and equilibrium line: (Ueq - Uop) = @(u(eq), current mt)
    diffLine = @(w, u, mt) s.abs2rel(u) - (s.ma./mt .* (s.abs2rel(w) - s.WR) + s.US);

    j = 0; % iteration count
    maxIterations = 50; tolerance = 1e-4; dmt = 0; % precision control
    while mt1-mt2 > s.mt*tolerance && j < maxIterations
        j = j + 1;
        % i is boolean value = whether there is an intersect, see detIntersect()
        i1 = detIntersect(s.aq, diffLine(s.aq, s.org, mt1), s.wR, s.wF);
        i2 = detIntersect(s.aq, diffLine(s.aq, s.org, mt2), s.wR, s.wF);

        % both values have intersection, we need to go higher
        if i1 && i2
            mt2 = mt1;
            mt1 = 10*mt1;
        % interval bisection, mt3 is the middle point
        elseif ~i1 && i2
            mt3 = 0.5*(mt1 + mt2);
            i3 = detIntersect(s.aq, diffLine(s.aq, s.org, mt3), s.wR, s.wF);
            if i3; mt2 = mt3; else mt1 = mt3; end
        else % this really shouldn't happen, but just in case...
            warning('Minimal mass of solvent cannot be calculated, no intersections have been found between operational and equilibrium line.');
        end
    end
    if j == maxIterations; warning('Minimal mass of solvent cannot be calculated, %d iterations exceeded!', j); end
    s.mSmin = 0.5*(mt1 + mt2) / (1 - s.uS);
    if s.mSmin >= s.mS; error('mS value specified or calculated from balance is lower than the minimal mS\nmS = %.4g\nmSmin = %.4g', s.mS, s.mSmin);end
end

% McCabe-Thiele method to get TP (number of theoretical plates)
function s = McCabeThiele(s)
    % operation line: U = @(W)
    opLine = @(W) s.ma/s.mt * (W - s.WR) + s.US;
    % initial value of W(alias x), U(alias y), TP and vectors for history of steps
    x = s.WF;
    y = s.UE;
    s.TP = 0;
    s.xV = []; s.yV = [];
    eqW = s.abs2rel(s.aq);
    eqU = s.abs2rel(s.org);

    % error because the first TP cannot be calculated - y is above the highest point of equilibrium curve and would require extrapolation (which yields disastrous effects in this case)
    % s.uE > max(s.aq) is equivalent to y > max(eqU) (which is *actually* what we need to check against)
    if s.uE > max(s.org); error('Cannot calculate theoretical plates, as uE exceeds experimental data range\nuE =  %.4g\nmax u(org) =  %.4g', s.uE, max(s.org)); end

    while x > s.WR
        % save the previous values
        s.xV = [s.xV x];
        s.yV = [s.yV y];

        % travel horizontally from operation line to equlibrium line: x = equilibrium(y)
        i = find(eqU > y, 1, 'first'); % i--
        x = eqW(i) + (y - eqU(i)) * (eqW(i-1) - eqW(i)) / (eqU(i-1) - eqU(i));

        % travel vertically from equlibrium line to operation line
        y = opLine(x);

        % increment TP, check for infinite loop
        s.TP = s.TP + 1;
        maxTP = 500; if s.TP == maxTP; error('Calculation of theoretical plates interrupted - reached %i and still not finished.', maxTP); end
    end

    s.xV = [s.xV x];
    s.yV = [s.yV y];

    % calculate the fraction of last TP by interpolation
    if x < s.WR; s.TP = s.TP - 1 + (s.WR - s.xV(end-1)) / (x - s.xV(end-1)); end
end

% detect intersect - returns true if the given vector U has a value that is lower than zero in such position that w is within bounds
function res = detIntersect(w, U, minVal, maxVal)
    res = false;
    for i = 1:length(w)
        if w(i) > minVal && w(i) < maxVal && U(i) < 0
            res = true; return;
        end
    end
end
