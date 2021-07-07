clear;
close all;

r = linspace(0, 5, 100);
rd = 0.5;
fatt = 1./((r/rd + 1).^2);

rd = 1;
fatt1 = 1./((r/rd + 1).^2);

figure();
plot(r, fatt, 'b-'); hold on;
plot(r, fatt1, 'g-');
plot(r, 1./(r.^2), 'r-');

grid on;
xlabel('Distance From Source (m)')
ylabel('Attenuation Factor (no units)')
ylim([0,1]);

legend('r_d = 0.5', 'r_d = 1.0', 'inverse square law');
