plot(input.diameter_range*1e6,phases/pi);
xlabel('Diameter [um]');
ylabel('Phase [rad]');
title("Wavelength " + num2str(input.lambda0) + ", period " + num2str(input.period_range) + ", thickness " + num2str(input.patterned_layer_t_range));
ylim([0, 2.05]);