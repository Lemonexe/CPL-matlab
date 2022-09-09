%s = DESIGNKARR(s, doDtCalc)
%       This function is a dependency for design.m and calculates geometry and operating variables for RECIPROCATING DISC COLUMNS (so called Karr columns)
%       State struct 's' is extended with newly calculated properties and returned
%       If 'doDtCalc' then only VCE and Dt is calculated and rest of calculation is skipped, otherwise, 's.Dt' is assumed and calculation proceeds
%       Calculation method is described in:
%           Holdup and Flooding in Reciprocating Plate Extraction Columns (1984) By H. M. I. Baird and Z. J. Shen
%               Psi, VCE, d32, phi
%           Computational Strategies for Sizing Liquid-Liquid Extractors (1999) by S. Hartland and A. Kumar
%               betaD, betaC
%           Axial Dispersion in a Reciprocating Plate Extraction Column (1976) by M. H. I. Baird and S. D. Kim
%               DAC

function s = designKarr(s, doDtCalc)
    addpath('lib'); % dependencies

    verb = s.verb >= 2 && ~doDtCalc; % whether console output will be emitted

    verb && headline('RECIPROCATING DISC COLUMN CALCULATIONS');

    g = 9.80665; % standard gravity [m/s2]
    drho = abs(s.rhoC-s.rhoD); % density difference [kg/m3]
    L = s.R;
    Af = s.stroke * s.N; % average linear velocity of plates [m/s]
    C0 = 0.6; % discharge coefficient through the perforations
    % C0 value taken from Flooding and Axial Dispersion in Reciprocating Plate Extraction Columns (1979) by M. M. Hafez, M. H. I. Baird and I. Nirdosh
    K = 30; % parameter of droplet circulation. 15 for perfectly circulating, 30 for rigid
    % 30 is more conservative - lower VCE

    phif = (sqrt(9*L^2 + 54*L + 1) - 7*L - 1) / (10 - 10*L); % retention on flooding (for flooding velocity, is not iterative)

    % ITERATIVE CALCULATIONS (for phi)
    % usually only two iterations are necessary
    phi = 0.05; % phi estimate
    phiprev = 1; % previous phi value
    i = 0; % iteration count
    maxIterations = 100; tolerance = 1e-4; % boundaries for while cycle
    while abs(phi/phiprev - 1) > tolerance && i < maxIterations
        i = i+1;
        phiprev = phi;

        rhomean = s.rhoC*(1-phi) + s.rhoD*phi; % mean density of dispersion [kg/m3]

        % this block is for VCE calculation, but also calculates Psi and d32, which is necessary for subsequent calculations
        Psi = 2*pi*pi/3 * ((1-s.aperture^2) / (s.H * C0^2 * s.aperture^2)) * rhomean * Af^3; % power dissipation per unit volume of dispersion [W/m3]
        Psim = Psi/rhomean; % Psi per mass [W/kg]
        d32 = 0.36 * s.tension^0.6 / (Psi^0.4 * rhomean^0.2); % drop size [m]
        uch = d32 * (g^2*drho^2/s.rhoC/s.muC)^(1/3) * K^(-2/3); % some kind of characteristic velocity [m/s]
        s.VCE = uch * (1-phif)^2 / (phif^(1/3)*(1-L) + L*phif^(-2/3)); % flooding velocity (cont) [m/s]

        % calculate s.Dt using correlations
        if doDtCalc
            VC = s.VCE * s.flowSafety; % superficial velocity of cont [m/s]
            VD = s.R * VC; % superficial velocity of disp [m/s]
            s.A = s.QC / VC; % column cross-sectional area [m2]
            s.Dt = sqrt(4*s.A/pi); % column diameter [m]
        % or use a fixed s.Dt value
        else
            s.A = pi * s.Dt^2 / 4;
            VC = s.QC / s.A;
            VD = s.R * VC;
        end

        callback = @(phi) VC - (1-phi)^2 * uch / ( phi^(1/3)*(1-L) + L*phi^(-2/3) );
        phi = fzero(callback, [1e-5, 0.3]);
        if isnan(phi); error('Couldn''t calculate phi, Dt is probably too big.'); end
    end

    % verb && fprintf(1, 'VCE = %.2e m/s (flooding velocity)\nVC =  %.2e m/s (velocity of cont phase at safety factor %.0f%%)\n', s.VCE, VC, s.flowSafety*100);
    % verb && fprintf(1, 'VD =  %.2e m/s (velocity of disp phase)\n', VD);
    % verb && fprintf(1, 'A = %.3f m2 (column cross-sectional area)\nDt = %.3f m (column diameter)\n\n', s.A, s.Dt);
    if doDtCalc; return; end % don't calculate anything else if the goal was to calculate Dt

    s.P = Psi * pi / 4 * s.Dt^2 * s.H; % power per disc (Psi * compartment volume)
    verb && fprintf(1, 'calculation finished in %i iterations at delta phi = %.3f%%\n', i, 100*abs(phi/phiprev - 1));
    verb && fprintf(1, 'P = %.2f W (power per disc)\nPsi = %.3f W/kg (gravimetric agitation power)\n\n', s.P, Psim);

    s.phi = phi;
    s.aint = 6*phi/d32; % interfacial surface [m2/m3]

    verb && fprintf(1, 'phi = %.3f (retention of disp in cont)\n\n', phi);
    verb && fprintf(1, 'd32 = %.2e m (mean drop diameter)\n', d32);
    verb && fprintf(1, 'a = %.2e m-1 (interfacial surface)\n\n', s.aint);

    % VslCorr = uch * (1-phi) / phi^(1/3); % direct correlation of slip velocity - is not needed, because it's actually equal to standard slip velocity formula:
    Vsl = VC/(1-phi) + VD/phi; % calculation of slip velocity from correlated retention
    verb && fprintf(1, 'Vsl = %.2e m/s\n\n', Vsl);

    % mass transfer coefs
    ReC = s.rhoC * d32 * Vsl / s.muC; % Reynolds in cont (this time it's Reynolds for mass transfer)
    verb && fprintf(1, 'ReC = %.1f (multiparticle Reynolds number)\n    should be between 10 and 500: %s\n', ReC, bool2str(ReC > 10 && ReC < 500));
    ScD = s.muD / s.rhoD / s.DD; % Schmidt in disp
    ScC = s.muC / s.rhoC / s.DC; % Schmidt in cont
    ShC_stag = 2.43 + ScC^(1/3)*(0.775*sqrt(ReC) + 0.0103*ReC); % constant
    ShC_circ = 50 + 2*sqrt(ReC * ScC / pi); % constant

    C1 = 2.44; %constant specific for column type
    kappa = s.muD/s.muC; % viscosity ratio
    aux1 = 1 + C1*( Psim/g * (s.rhoC/g/s.tension)^(1/4) )^(1/3);
    aux2 = 5.26e-2 * ReC^(1/3 + 6.59e-2*ReC^0.25) * (ScC * Vsl * s.muC / s.tension)^(1/3) * 1/(1+kappa^1.1);
    aux3 = 17.7 + (3.19e-3*(ReC*ScD^(1/3))^1.7) / (1+1.43e-2*(ReC*ScD^(1/3))^0.7) * (s.rhoD/s.rhoC)^(2/3) * 1/(1+kappa^(2/3));

    ShC = (aux2*aux1*ShC_circ + ShC_stag) * (1-phi) / (1+aux2*aux1);

    ShD = aux3*aux1;

    s.betaD = ShD * s.DD / d32;
    s.betaC = ShC * s.DC / d32;
    verb && fprintf(1, 'ShD = %.2e (Sherwood in disp)\nShC = %.2e (Sherwood in cont)\n', ShD, ShC);
    verb && fprintf(1, 'betaD = %.2e m/s (mass transfer coef in disp)\nbetaC = %.2e m/s (mass transfer coef in cont)\n\n', s.betaD, s.betaC);

    % axial dispersion - empirical formula only available for continuous, therefore only Peclet number used for dispersed
    s.DAC = 5.62 * (s.stroke*50)^1.41 * s.N^0.73 * (s.H*100)^-0.88 * 1e-4;
    s.DAD = s.DAC*VD/VC; % PeD = PeC
    verb && fprintf(1, 'DAD = %.2e m2/s (axial dispersion coef in disp)\nDAC = %.2e m2/s (axial dispersion coef in cont)\n\n', s.DAD, s.DAC);
end
