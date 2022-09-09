%s = DESIGNROTDISC(s, doDtCalc)
%       This function is a dependency for design.m and calculates geometry and operating variables for ROTATING DISC or KÜHNI EXTRACTORS
%       State struct 's' is extended with newly calculated properties and returned
%       If 'doDtCalc' then only VCE and Dt is calculated and rest of calculation is skipped, otherwise, 's.Dt' is assumed and calculation proceeds
%       Calculation method is described in "Design of Liquid-Liquid Extractors" (1997) by S. Hartland and A. Kumar

function s = designRotDisc(s, doDtCalc)
    addpath('lib'); % dependencies

    verb = s.verb >= 2 && ~doDtCalc; % whether console output will be emitted

    g = 9.80665; % standard gravity [m/s2]
    drho = abs(s.rhoC-s.rhoD); % density difference [kg/m3]
    C2D = s.dispPhase == 2; % whether mass transfer takes place from continuous to dispersed

    % struct object with constants that depend on column type
    c = struct();
    if strcmpi(s.type, 'rotDisc')
        verb && headline('ROTATING DISC COLUMN CALCULATIONS');
        c.e = 0.49; % fractional free cross-sectional area = (DS/DC)^2, is by definiton equal to 0.7^2
        c.CR = c.e; % constant for flooding velocity
        c.DR = @(DC) 0.6 * DC;
        c.DS = @(DC) 0.7 * DC;
        c.H  = @(DC) 0.13 * DC^0.67;
        c.NP = @(ReR) 109.36/ReR + 0.74 * ((1e3 + 1.2*ReR^0.72) / (1e3 + 3.2*ReR^0.72))^3.3;
        % various unnamed numeric constants
        c.c = [1, 1.29, 2.54, 0.97, 0.64, -0.45, -1.12, 0.19, 1.06, 0.95, 0.67, 0.69, 7.13, -0.65, 0.14, -0.26, -0.10];
        c.C3 = @(DR,H) (DR / H)^0.62;
        c.betaD = @(Vsl, ScD) 0.023 * Vsl / sqrt(ScD);
    elseif strcmpi(s.type, 'kuhni')
        verb && headline('KÜHNI COLUMN CALCULATIONS');
        c.e = 0.3; % ranges from 0.2 to 0.4
        c.CR = 0.25 + 0.75*c.e;
        c.DR = @(DC) 5/12 * DC; % ranges from 0.33 to 0.5
        c.DS = @(DC) NaN; % not applicable for a Kühni column
        c.H  = @(DC) 0.25 * DC^0.6; % ranges from 0.2 to 0.3
            % c.DR = @(DC) 0.5 * DC; c.DS = @(DC) NaN; c.H  = @(DC) 0.3 * DC^0.6; % AN ALTERNATE SETTING FOR VALIDATION!
        c.NP = @(ReR) 1.08 + 10.94/sqrt(ReR) + 257.37/(ReR^1.5);
        c.c = [1, 3.04, 1.60, 3.40e-2, 0.45, -0.63, -0.38, 0.0267, 1, 0.56, 0.77, 0.64, 20.7, -0.34, 0, -0.77, 0];
        c.C3 = @(DR,H) 2.27;
        c.betaD = @(Vsl, ScD) 0.00375 * Vsl / (1 + s.muD/s.muC);
    end


    % ITERATIVE CALCULATIONS (for Dt)
    % only if Dt is calculated, otherwise the cycle is broken
    phie = 0;
    VC = 4e-3; % initial estimate of continuous superficial velocity [m/s]
    DCprev = 0; DC = 1e10;
    i = 0; % iteration count
    maxIterations = 100; tolerance = 1e-4; % boundaries for while cycle
    while abs(DC/DCprev - 1) > tolerance && i < maxIterations
        i = i+1;

        % calculate s.Dt and column geometry using correlations
        if doDtCalc
            s.A = s.QC / VC; % column cross-sectional area [m2]
            DCprev = DC; % save the old value of DC before calculating new one (while cycle condition)
            DC = sqrt(4*s.A/pi); % column diameter [m]
        % or use a fixed s.Dt value
        else
            DC = s.Dt;
            s.A = pi * s.Dt^2 / 4;
            VC = s.QC / s.A;
            VD = s.R * VC;
        end

        % other geometry
        DR = c.DR(DC); % rotor diameter [m]
        DS = c.DS(DC); % clearance diameter [m]
        H = c.H(DC); % compartment height [m]

        % power characteristics
        ReR = s.N * DR*DR * s.rhoC / s.muC; % Reynolds of rotating disc
        NP = c.NP(ReR); % power number
        P = NP * s.N^3 * DR^5 * s.rhoC; % power per disc [W]
        Psi = 4 * P / (pi * DC*DC * H * s.rhoC); % gravimetric agitation power [W/kg]

        % drop size (for transfer coefs)
        if C2D; C1 = c.c(1); else C1 = c.c(2); end
        C2 = c.c(3); C3 = c.c(4); n1 = c.c(5); n2 = c.c(6); n3 = c.c(7);
        aux1 = C2 * sqrt(s.tension / (drho * g * H^2));
        aux2 = C3 * (Psi/g * (s.rhoC/g/s.tension)^0.25)^n2 * (H * sqrt(s.rhoC*g/s.tension))^n3;
        d32 = H * C1 * c.e^n1 / (1/aux1 + 1/aux2);

        % retention on flooding (for flooding velocity)
        C1 = 24 * s.muC / (0.53 * d32 * s.rhoC);
        C2 = 4 * d32 * g * drho / (1.59 * s.rhoC);
        C3 = @(phie) sqrt(C1^2 + 4*C2*(1-phie)/(1+4.56*phie^0.73));

        aux1 = @(phie) 2 * C2 * phie * (1-phie) / (C3(phie) * (1+4.56*phie^0.73)^2);
        aux2 = @(phie) 1 + 4.56*phie^0.73 + 3.33*phie^-0.27*(1+phie);
        aux3 = @(phie) (C3(phie) - C1) * phie * (1-phie) * (s.R-1);

        callback = @(phie) ((C3(phie)-C1)*(1-2*phie) - aux1(phie)*aux2(phie)) * (phie + s.R*(1-phie)) + aux3(phie);
        phie = fzero(callback, 0.2); % solve equation with initial estimate of retention = 0.2

        % if Dt was specified, do not calculate flooding and break the cycle to continue with currently obtained d32 and Psi etc.
        if ~doDtCalc; break; end

        % flooding velocity (for determination of optimal superficial velocities)
        s.VCE = c.CR * (-C1 + C3(phie)) * phie * (1-phie) / 2 / (phie + s.R*(1-phie)); % flooding velocity (cont) [m/s]
        VC = s.VCE * s.flowSafety; % superficial velocity of cont [m/s]
        VD = VC * s.R; % superficial velocity of disp [m/s]
    end
    s.Dt = DC; s.H = H; s.P = P;
    if doDtCalc; return; end % don't calculate anything else if the goal was to calculate Dt

    % verb && fprintf(1, 'calculation finished in %i iterations at deltaDC = %.3f%%\n', i, 100*abs(DC/DCprev - 1));
    % verb && fprintf(1, 'A = %.3f m2 (column cross-sectional area)\nDC = %.3f m (column diameter)\nDR = %.3f m (disc diameter)\nDS = %.3f m (clearance diameter)\nH = %.3f m (compartment height)\n\n', s.A, DC, DR, DS, H);
    verb && fprintf(1, 'P = %.2f W (power per disc)\nPsi = %.3f W/kg (gravimetric agitation power)\nN = %.2f Hz (frequency)\nReR = %i (Reynolds number of mixing)\n', P, Psi, s.N, round(ReR));
    verb && fprintf(1, 'd32 = %.2e m (mean drop diameter)\n    should be between 1.0e-3 and 1.5e-3: %s\n\n', d32, bool2str(d32 > 1e-3 && d32 < 1.5e-3));
    % verb && fprintf(1, 'VCE = %.2e m/s (flooding velocity)\nVC =  %.2e m/s (velocity of cont phase at safety factor %.0f%%)\n', s.VCE, VC, s.flowSafety*100);
    % verb && fprintf(1, 'VD =  %.2e m/s (velocity of disp phase)\n\n', VD);


    % retention (for slip velocity & transfer coefs)
    C1 = c.c(8);
    if C2D; C2 = c.c(9); else C2 = c.c(10); end
    C3 = c.C3(DR, H);
    n1 = c.c(11); n2 = c.c(12); n3 = c.c(13); n4 = c.c(14); n5  = c.c(15); n6 = c.c(16); n7 = c.c(17);
    aux0 = (s.rhoC / g / s.tension)^0.25;
    aux1 = (C1 + (Psi/g*aux0)^n1) * (VD * aux0)^n2 * exp(n3 * VC * aux0);
    muW = 1e-3; % standard viscosity (water @20°C) [Pa*s]
    aux2 = C2 * (drho/s.rhoC)^n4 * (s.muD/muW)^n5 * C3 * c.e^n6 * (H * sqrt(s.rhoC*g/s.tension))^n7;

    phi = aux1 * aux2;
    if phi >= phie && phie > 0; warning('retention of disp in cont (%.2f) has exceeded flooding retention (%.2f). Use lower frequency.', phi, phie); end
    s.phi = phi;
    s.aint = 6*phi/d32; % interfacial surface [m2/m3]
    Vsl = VC/(1-phi) + VD/phi; % calculation of slip velocity from correlated retention
    verb && fprintf(1, 'phi = %.3f (retention of disp in cont)\nVsl = %.2e m/s (slip velocity)\n\n', phi, Vsl);
    verb && fprintf(1, 'a = %.2e m-1 (interfacial surface)\n\n', s.aint);

    % axial dispersion
    aux0 = s.N*DR/VC;
    if strcmpi(s.type, 'rotDisc')
        aux1 = (1.26e-2*aux0 + 13.38/(3.18+aux0)) * (VC*DR*s.rhoC/s.muC)^-0.08 * (DC/DR)^0.16 * (DC/H)^0.10 * c.e;
        s.DAC = VC * H * (0.42 + 0.29*s.R + aux1);
        aux2 = 9.37*(s.N*DR*phi/VD)^-0.48 * (DR^2*drho*g/s.tension)^-0.64 *(DC/DR)^0.70 * phi^-0.9;
        s.DAD = VD * DR * (0.30 * (VC/(1-phi)+VD/phi) / (VD/phi) + aux2);
    elseif strcmpi(s.type, 'kuhni')
        s.DAC = VC * H * (0.14 + 0.046*aux0*phi*c.e);
        s.DAD = VD * H * (22.62 * c.e^0.43 * (s.N*DR^2*phi/VD/DC)^-0.26 * (d32*VD*s.rhoD/phi/s.muD)^-0.57);
    end
    verb && fprintf(1, 'DAD = %.2e m2/s (axial dispersion coef in disp)\nDAC = %.2e m2/s (axial dispersion coef in cont)\n\n', s.DAD, s.DAC);

    % mass transfer coefs
    ReC = s.rhoC * d32 * Vsl / s.muC; % Reynolds in cont (this time it's Reynolds for mass transfer)
    verb && fprintf(1, 'ReC = %.1f (multiparticle Reynolds number)\n    should be between 10 and 500: %s\n', ReC, bool2str(ReC > 10 && ReC < 500));
    ScD = s.muD / s.rhoD / s.DD; % Schmidt in disp
    ScC = s.muC / s.rhoC / s.DC; % Schmidt in cont
    ShC_stag = 2.43 + ScC^(1/3)*(0.775*sqrt(ReC) + 0.0103*ReC); % constant
    ShC_circ = 2*sqrt(ReC * ScC / pi); % constant
    ShC = ShC_stag + (ShC_circ - ShC_stag) * (1 - exp(-4.18e-3*(ReC*ScC)^0.42)); % Sherwood in cont
    s.betaC = ShC * s.DC / d32;
    s.betaD = c.betaD(Vsl, ScD);
    verb && fprintf(1, 'ShC = %.2e (Sherwood in cont)\n', ShC);
    verb && fprintf(1, 'betaD = %.2e m/s (mass transfer coef in disp)\nbetaC = %.2e m/s (mass transfer coef in cont)\n\n', s.betaD, s.betaC);
end
