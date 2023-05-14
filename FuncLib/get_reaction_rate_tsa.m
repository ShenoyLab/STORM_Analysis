%% This function generate the code extracting the statistical propertity of
% chemical reaction rate
% input: data --- old
% output: data --- old
function data = get_reaction_rate_tsa(data)

num_of_data = length(data);

% get number of treatments
names = string([]);
for cell_idx = 1:num_of_data
%     names = [names;extractBefore(data{cell_idx,1}.name,'_')];
    names = [names;data{cell_idx,1}.name(1:find(data{cell_idx,1}.name == '_', 1, 'last')-1)];
end
num_of_treatments = length(unique(names));
sample_info = zeros(num_of_treatments,1);
temp = unique(names);
for i = 1:num_of_treatments
    sample_info(i) = length(names(names==temp(i)));
end

for cell_idx = 1:num_of_data
    monitor = figure('name',data{cell_idx,1}.name);
    nbins = 70;
    radius = data{cell_idx,1}.hetero_radius; R = radius/1000;
    D = 1;
    
    spacing = data{cell_idx,1}.spacing; d = spacing/1000;
    
    subplot(2,2,1)% Radius
    [p,x_bins] = hist(R,nbins);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',2),'b','linewidth',2); %PDF
    title('PDF for Radius');xlabel('Radius /\mum');
    
    subplot(2,2,2)% Spacing
    [p,x_bins] = hist(d,nbins);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',2),'b','linewidth',2); %PDF
    title('PDF for Spacing');xlabel('Spacing /\mum');
    
    subplot(2,2,3)% G_ac Fitting
    g_ac = D./(d.^2);
    [p_ac,x_bins_ac] = hist(g_ac,nbins);
    p_ac = p_ac/sum(p_ac)/(x_bins_ac(2)-x_bins_ac(1));
    [xData, yData] = prepareCurveData( x_bins_ac, p_ac );

    % Set up fittype and options.
    ft = fittype( 'exp(-(log(x)-mu)^2/(2*sigma^2))/(x*sigma*sqrt(2*pi))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [3 0.2];

    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    data{cell_idx,1}.g_ac_coeffvals = coeffvalues(fitresult{1});
    % Plot fit with data.
    
    h = plot( fitresult{1}, xData, yData );
    legend( h, '\Gamma_{ac}', 'Fitting pdf', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( '\Gamma_{ac}', 'Interpreter', 'none' );
    ylabel( 'pdf', 'Interpreter', 'none' );
    grid on
    
    subplot(2,2,4)% G_me Fitting
%     g_me = (R.^2).*(g_ac.^2)./(D-(R.^2).*(g_ac));
    g_me = g_ac./((d./R).^2-1);
    g_me = g_me(g_me>0);
    g_me = g_me(g_me<=20);

    [p_me,x_bins_me] = hist(g_me,15);
    p_me = p_me/sum(p_me)/(x_bins_me(2)-x_bins_me(1));
    [xData, yData] = prepareCurveData( x_bins_me, p_me );

    % Set up fittype and options.
    ft = fittype( 'exp(-(log(x)-mu)^2/(2*sigma^2))/(x*sigma*sqrt(2*pi))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [-5 3];

    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    data{cell_idx,1}.g_me_coeffvals = coeffvalues(fitresult{1});
    % Plot fit with data.
    
    h = plot( fitresult{1}, xData, yData );
    legend( h, '\Gamma_{me}', 'Fitting pdf', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( '\Gamma_{me}', 'Interpreter', 'none' );
    ylabel( 'pdf', 'Interpreter', 'none' );
    grid on
    pause(1)
    close(monitor)
    
end
% ----------------------------------------------------------------------
% ----------------------------- SUMMARY --------------------------------
% ----------------------------------------------------------------------

figure('name','summary')
subplot(321)
for cell_idx = 1:num_of_data
    if cell_idx<=6
        cmap = 'r';
    else
        cmap = 'b';
    end
    spacing = data{cell_idx,1}.spacing; d = spacing/1000;
    [p,x_bins] = hist(d,nbins);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',2),cmap,'linewidth',2); hold on
    title('PDF for Spacing');xlabel('Spacing /\mum');
end

subplot(322)
for cell_idx = 1:num_of_data
    if cell_idx<=6
        cmap = 'r';
    else
        cmap = 'b';
    end
    radius = data{cell_idx,1}.hetero_radius; R = radius/1000;
    [p,x_bins] = hist(R,nbins);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',2),cmap,'linewidth',2); hold on
    title('PDF for Radius');xlabel('Radius /\mum');
    xlim([0,0.3]);
end

subplot(323)
g_ac_mean = [];
for cell_idx = 1:num_of_data
    if cell_idx<=sample_info(1)
        cmap = 'r';
    elseif cell_idx <= sample_info(1)+sample_info(2)
        cmap = 'b';
    else
        cmap = 'g';
    end
    
    mu = data{cell_idx,1}.g_ac_coeffvals(1);
    sigma = data{cell_idx,1}.g_ac_coeffvals(2);
    g_ac_mean = [g_ac_mean;exp(mu+sigma^2/2)];
    x_hat = linspace(0,80,1000);
    p_hat = exp(-(log(x_hat)-mu).^2/(2.*sigma.^2))./(x_hat.*sigma.*sqrt(2.*pi));
    plot(x_hat, p_hat,cmap,'linewidth',1.5);hold on
    title('PDF for \Gamma_{ac}');xlabel('\Gamma_{ac} /\mum');
end
subplot(324)
g_me_mean = [];
for cell_idx = 1:num_of_data
    if cell_idx<=6
        cmap = 'r';
    else
        cmap = 'b';
    end
    
    mu = data{cell_idx,1}.g_me_coeffvals(1);
    sigma = data{cell_idx,1}.g_me_coeffvals(2);
    g_me_mean = [g_me_mean;exp(mu+sigma^2/2)];
    x_hat = linspace(0,3,10000);
    p_hat = exp(-(log(x_hat)-mu).^2/(2.*sigma.^2))./(x_hat.*sigma.*sqrt(2.*pi));
    plot(x_hat, p_hat,cmap,'linewidth',1.5);hold on
%     xlim([0,10]);
    ylim([0,2]);
    title('PDF for \Gamma_{me}');xlabel('\Gamma_{me} /\mum');
end
subplot(325)
ac_mean = [mean(g_ac_mean(1:6)),mean(g_ac_mean(7:11))];
ac_std = [std(g_ac_mean(1:6)),std(g_ac_mean(7:11))];
bar(ac_mean); hold on
errorbar(1:2,ac_mean,ac_std,ac_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
set(gca,'xticklabel',{'Control';'TSA'})

subplot(326)
me_mean = [mean(g_me_mean(1:6)),mean(g_me_mean(7:11))];
me_std = [std(g_me_mean(1:6)),std(g_me_mean(7:11))];
bar(1:2,me_mean); hold on
errorbar(1:2,me_mean,me_std,me_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
set(gca,'xticklabel',{'Control';'TSA'})
end