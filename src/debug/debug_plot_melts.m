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

    x = frac.mesh.xc;
    plot(x, frac.betaeq, LineWidth=2)
    hold on;
    plot(x, frac.beta, LineWidth=2)
    plot(x, frac.crystallization.betaeqList{1}, LineWidth=2)
    plot(x, frac.crystallization.betaList{1}, LineWidth=2)
    plot(x, frac.crystallization.betaeqList{2}, LineWidth=2)
    plot(x, frac.crystallization.betaList{2}, LineWidth=2)
    legend("$\beta_{eq}$", "$\beta$", "$\beta_{eq, plag}$", "$\beta_{plag}$", "$\beta_{eq, mafic}$", "$\beta_{mafic}$")
end