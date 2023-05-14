function coeff = logNfit(gamma,nbins,guess)
    
    [p,x_bins] = hist(gamma,nbins);
    p = p/sum(p)/(x_bins(2)-x_bins(1));
    [xData, yData] = prepareCurveData( x_bins, p);

    % Set up fittype and options.
    ft = fittype( 'exp(-(log(x)-mu)^2/(2*sigma^2))/(x*sigma*sqrt(2*pi))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = guess;
    
    % Fit model to data.
    [fitresult{1}, log] = fit( xData, yData, ft, opts );
    coeff = coeffvalues(fitresult{1});
    if log.rsquare < 0.8
        fprintf('Bad Fit! \n')
    end
    % Plot fit with data.
    monitor = figure();
    h = plot( fitresult{1}, xData, yData );
    legend( h, '\Gamma', 'Fitting pdf', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( '\Gamma_{ac}', 'Interpreter', 'none' );
    ylabel( 'pdf', 'Interpreter', 'none' );
    grid on
    pause(1)
    close(monitor)
end