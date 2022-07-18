function [S,cb] = plotlinescanviewplane(Y, Z, C, maxC, T_S_2_LS, T_F_2_LS, viewPlaneFrustrum, lsFigYlims, lsFigZlims, cm)
% Plot the data present on the view-plane of a line-scan camera. View-plane
% is assuemd to lie on the YZ plane. Can be used to plot irradiance,
% variance, absolute difference.
arguments
    Y
    Z
    C
    maxC
    T_S_2_LS
    T_F_2_LS
    viewPlaneFrustrum
    lsFigYlims
    lsFigZlims
    cm = hot(1000);
end

%plot surface
S = surf(Y, Z, C, 'EdgeColor', 'none'); hold on;
colormap(cm);
cb = colorbar('Ticks', 0:maxC/5:maxC);
caxis([0, maxC]);
xlim(lsFigYlims);
ylim(lsFigZlims);

%line-scan camera position (origin) (blue)
scatter3(0, 0, (maxC*10), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); 

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxC*10), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame (red)
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxC*10), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
yPoly = viewPlaneFrustrum(1,:);
zPoly = viewPlaneFrustrum(2,:);
plot3(yPoly, zPoly, (maxC*10)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

xlabel('y (m)');
ylabel('z (m)');
grid off;
view(2);
hold off;

end