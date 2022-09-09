%s = DESIGNPACKED(s, doDtCalc)
%       This function is a dependency for design.m and calculates geometry and operating variables for PACKED COLUMNS
%       State struct 's' is extended with newly calculated properties and returned
%       If 'doDtCalc' then only VCE and Dt is calculated and rest of calculation is skipped, otherwise, 's.Dt' is assumed and calculation proceeds
%       Calculation method is described in:
%           "Liquid-Liquid and Solid-Liquid Extractors" (2016) by Jean-Paul Duroudier
%               phi, d32, Vsl, betaD, betaC, DAC, DAD, Dncrit
%           "Design of Liquid-Liquid Extractors" (1997) by S. Hartland and A. Kumar
%               VCE

function s = designPacked(s, doDtCalc)
    addpath('lib'); % dependencies

    verb = s.verb >= 2 && ~doDtCalc; % whether console output will be emitted

    verb && headline('PACKED COLUMN CALCULATIONS');

    g = 9.80665; % standard gravity [m/s2]
    drho = abs(s.rhoC-s.rhoD); % density difference [kg/m3]
    C2D = s.dispPhase == 2; % whether mass transfer takes place from continuous to dispersed

    % unpack packing struct
    packingType = s.packing.type;
    Dn = s.packing.Dn;
    at = s.packing.at;
    eps = s.packing.eps;

    % calculate only s.Dt using correlations and return
    if doDtCalc
        % flooding velocity (for determination of optimal superficial velocities)
        auxStruct = struct('Raschig', 0.28, 'Intalox', 0.28, 'Berl', 0.37, 'Lessing', 0.37, 'Pall', 0.20, 'Bialecki', 0.25, 'sphere', 0.30);
        if ~isfield(auxStruct, packingType);
            warning(['Packing type ''' packingType ''' not found, defaulting to Raschig']);
            packingType = 'Raschig';
        end
        C1 = auxStruct.(packingType);
        aux1 = sqrt(g/at) / (1 + sqrt(s.R))^2 * C1 * eps^1.54 * (drho/s.rhoD)^0.41;
        aux2 = (1/at*(drho^2*g/s.muC^2)^(1/3))^0.3;
        aux3 = (s.muC*sqrt(at/drho/s.tension))^0.15;
        s.VCE = aux1*aux2*aux3; % flooding velocity (cont) [m/s]
        VC = s.VCE * s.flowSafety; % superficial velocity of cont [m/s]
        VD = s.R * VC; % superficial velocity of disp [m/s]
        s.A = s.QC / VC; % column cross-sectional area [m2]
        s.Dt = sqrt(4*s.A/pi); % column diameter [m]
        % verb && fprintf(1, 'VCE = %.2e m/s (flooding velocity)\nVC =  %.2e m/s (velocity of cont phase at safety factor %.0f%%)\n', s.VCE, VC, s.flowSafety*100);
        % verb && fprintf(1, 'VD =  %.2e m/s (velocity of disp phase)\n', VD);
        % verb && fprintf(1, 'A = %.3f m2 (column cross-sectional area)\nDt = %.3f m (column diameter)\n\n', s.A, s.Dt);
        return
    end

    % or use a fixed s.Dt value and proceed
    s.A = pi * s.Dt^2 / 4;
    VC = s.QC / s.A;
    VD = s.R * VC;

    % proper size of packing
    Dncrit = 2.42 * sqrt(s.tension / drho / g); % critical size of packing
    Dnmax = s.Dt / 10;
    verb && fprintf(1, 'Size of packing = %.2e m\nmin = %.2e m (critical)\nmax = %.2e m (Dt/10)\nis within limits: %s\n\n', Dn, Dncrit, Dnmax, bool2str(Dn >= Dncrit && Dn <= Dnmax));

    % retention (for slip velocity & transfer coefs)
    if C2D; C1 = 6.16; else C1 = 3.76; end
    aux1 = C1 * eps^-1.11 / sqrt(drho / s.rhoC);
    aux2 = (1/at * (s.rhoC^2*g/s.muC^2)^(1/3))^-0.72;
    aux3 = (s.muD/s.muC)^0.1 * (VD * (s.rhoC/g/s.muC)^(1/3))^1.03;
    aux4 = exp(0.95*VC*(s.rhoC/g/s.muC)^(1/3));
    phi = aux1*aux2*aux3*aux4;
    s.phi = phi;
    verb && fprintf(1, 'phi = %.3f (retention of disp in cont)\n\n', phi);

    % drop size (for transfer coefs)
    if C2D; C1 = 2.24; else C1 = 3.13; end
    d32 = sqrt(s.tension/g/drho) * C1 * (s.muC*g^0.25*drho^-0.25*s.tension^-0.75*s.rhoC/s.rhoD)^0.19;
    s.aint = 6*eps*phi/d32; % interfacial surface [m2/m3]
    verb && fprintf(1, 'd32 = %.2e m (mean drop diameter)\n', d32);
    verb && fprintf(1, 'a = %.2e m-1 (interfacial surface)\n\n', s.aint);

    % slip velocity (for transfer coefs)
    if C2D; C1 = 0.27; else C1 = 0.38; end
    aux1 = (g*s.muC/s.rhoC)^(1/3) * C1 * eps^-0.17;
    aux2 = (drho/s.rhoC)^0.41 * (1/at * (s.rhoC^2*g/s.muC^2)^(1/3))^0.59;
    aux3 = (s.muD/s.muC)^-0.1 * (1 - phi);
    VslCorr = aux1*aux2*aux3; % direct correlation of slip velocity
    VslCalc = VC/eps/(1-phi) + VD/eps/phi; % calculation of slip velocity from correlated retention
    VslMean = 0.5*VslCorr + 0.5*VslCalc;
    verb && fprintf(1, 'Vsl = %.2e m/s (slip velocity correlated)\nVsl = %.2e m/s (slip velocity calculated)\nVsl = %.2e m/s (slip velocity mean)\n\n', VslCorr, VslCalc, VslMean);

    % transfer coefs
    ReC = s.rhoC * d32 * VslMean / s.muC; % Reynolds in cont
    ScD = s.muD / s.rhoD / s.DD; % Schmidt in disp
    ScC = s.muC / s.rhoC / s.DC; % Schmidt in cont

    s.betaD = 0.023 * VslMean * ScD^-0.5;

    ShC = 0.698 * ScC^0.4 * sqrt(ReC) * (1-phi); % Sherwood in cont
    s.betaC = ShC * s.DC / d32; % mass transfer coef from Sherwood
    verb && fprintf(1, 'ShC = %.2e (Sherwood in cont)\n', ShC);
    verb && fprintf(1, 'betaD = %.2e m/s (mass transfer coef in disp)\nbetaC = %.2e m/s (mass transfer coef in cont)\n\n', s.betaD, s.betaC);

    % dispersion coefs
    PeD = 0.45*exp(-0.405 / s.R); % Peclet number
    s.DAD = VD * Dn / PeD; % axial dispersion coef for disp

    aux1 = s.R * sqrt(0.6*s.muC/s.rhoC/Dn/VC);
    PeC = 0.01268*aux1^-0.9043;
    ReC = VC*s.rhoC*Dn/s.muC/(1-eps); % Reynolds in cont
    PecLim = 0.1443*log(ReC) - 0.232;
    if PeC > PecLim; PeC = PecLim; end % Pec has to be smaller than PecLim, otherwise we will use PecLim instead of Pec
    s.DAC = VC*Dn*eps/PeC; % axial dispersion coef for cont
    verb && fprintf(1, 'DAD = %.2e m2/s (axial dispersion coef in disp)\nDAC = %.2e m2/s (axial dispersion coef in cont)\n\n', s.DAD, s.DAC);
end
