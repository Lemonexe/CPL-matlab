%EXTRACTORPLOT(v, sol, SA, t)
%       see the extractor.m function - it calculates and returns the raw data
%       this function serves to write the extractor() results into console and draw a plot of McCabe-Thiele procedure
%
%CALL EXAMPLE:
%       extractorPlot([100, 10, 0.01, 0.001, 0, NaN], 'TOL', 40, 60)

function extractorPlot(v, sol, SA, t)
    % call the main function to get raw data
    s = extractor(v, sol, SA, t);

    % text outputs
    varNames = {'mF', 'mS', 'wF', 'wR', 'uS', 'uE'}; %labels for the variables in 'v' vector
    fprintf(1, '\nMass balance:\n%s = %.4g\nmR = %.4g\nmE = %.4g\n', varNames{s.i}, s.v(s.i), s.mR, s.mE);
    fprintf(1, '\nminimal mS = %.4g\n', s.mSmin);
    fprintf(1, '\nTheoretical plates = %.2f\n', s.TP);

    % initiate plot, draw operational line
    figure();
    plot([s.WF s.WR], [s.UE s.US], ':xr');
    %GRAYSCALE plot([s.WF s.WR], [s.UE s.US], '--xk');
    hold on;
    xlabel('\itW\rm_{CPL, aq}');
    ylabel('\itW\rm_{CPL, org}');
    title(sprintf('%s, %d %% SA, %d °C', s.sol, s.SA, s.t));

    % plot equilibrium curve from from getEqFun
    plot(s.abs2rel(s.aq), s.abs2rel(s.org), '-k');
    %GRAYSCALE plot(s.abs2rel(s.aq), s.abs2rel(s.org), '-.k');
    axis([0 s.WF*1.2 0 s.UE*1.2]);

    % plot all theoretical plates
    xPlot = zeros(1, 2*6-1); yPlot = xPlot;
    for i = 1:(length(s.xV) - 1)
        xPlot([2*i-1,2*i]) = [s.xV(i), s.xV(i+1)];
        yPlot([2*i-1,2*i]) = [s.yV(i), s.yV(i)];
    end
    xPlot(end) = s.xV(end);
    yPlot(end) = s.yV(end);
    plot(xPlot, yPlot, 'g-');
    %GRAYSCALE plot(xPlot, yPlot, '-k');

    hold off;
end
