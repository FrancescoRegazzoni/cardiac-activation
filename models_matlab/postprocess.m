% Author: Francesco Regazzoni - MOX, Politecnico di Milano
% Email:  francesco.regazzoni@polimi.it
% Date:   2020

function postprocess(output, filename)

    figure('units','pixel','outerposition', [100, 100, 600, 400]);

    subplot(2,2,1);
    plot(output.t, output.Ca, 'linewidth', 2)
    ylim([0 1.2])
    ylabel('[Ca^{2+}]_i [\muM]')

    subplot(2,2,2);
    plot(output.t, output.SL, 'linewidth', 2)
    ylim([2 2.3])
    ylabel('SL [\mum]')

    subplot(2,2,3);
    plot(output.t, output.P, 'linewidth', 2)
    ylim([0 1])
    ylabel('permissivity [-]')

    subplot(2,2,4);
    plot(output.t, output.Ta, 'linewidth', 2)
    ylim([0 80])
    ylabel('active tension [kPa]')

    for i = 1:4
        subplot(2,2,i);
        xlabel('time [s]')
        xlim([-inf inf])
    end

    set(gcf, 'PaperPosition', [0 0 6 4]);
    set(gcf, 'PaperSize', [6 4]);
    saveas(gcf, filename)

end