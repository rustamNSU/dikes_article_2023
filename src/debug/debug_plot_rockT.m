function debug_plot(frac)
    stFontSize = 12;
    set(0, 'defaultTextInterpreter','latex')
    set(0, 'defaultTextFontName', 'CMU Serif')
    set(0, 'defaultAxesFontName', 'CMU Serif')
    set(0, 'defaultAxesFontSize', stFontSize)
    set(0, 'defaultAxesTickLabelInterpreter', 'latex')
    set(0, 'defaultLineLineWidth', 2)
    set(0, 'defaultLegendFontName', 'CMU Serif')
    set(0, 'defaultLegendFontSize', stFontSize)
    set(0, 'defaultLegendInterpreter', 'latex')

    xc = frac.mesh.xc;
    yc = frac.rockTemperature.yc;
    [X, Y] = meshgrid(xc, yc);
    Z = frac.rockTemperature.T' - 273.15;
    contourf(Y, X, Z, LevelList=linspace(min(Z, [], 'all'), max(Z, [], 'all'), 50), EdgeColor='none');
    colormap('turbo')
    h_cb = colorbar();
    ylabel(h_cb, "Temperature, C$^\circ$", FontSize=stFontSize, Interpreter='latex', Visible='on');
    set(h_cb, TickLabelInterpreter='latex', FontSize=stFontSize);
end