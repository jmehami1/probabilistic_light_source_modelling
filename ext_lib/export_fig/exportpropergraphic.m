function exportpropergraphic(fig, imgName, hwRatio, fontSize, picturewidth)
% Properly formats a MATLAB figure that will be exported to image format for
% research papers, report, and thesis
% INPUTS:
%       fig - figure which will be exported
%       imgName - fullpath, filename and filetype of exported figure image
%       hwRatio - height to width ratio of exported image (height/weight)
%       fontSize - font size for all text in figures. Note that some label
%           texts are not affected by this change E.g polarplots and colorbar titles
%       picturewidth - width of exported image

arguments
    fig
    imgName
    hwRatio = 1
    fontSize = 18
    picturewidth = 20
end

%should set all axes and tick colors to black
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});
fprintf('Exporting %s\n', imgName);

set(findall(fig, '-property', 'FontSize'),'FontSize',fontSize); % adjust fontsize to your document
set(findall(fig, '-property', 'FontName'),'FontName','times'); % adjust fontsize to your document
set(findall(fig,'-property','Interpreter'),'Interpreter','latex');
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig,'Units','centimeters','Position',[3 3 picturewidth hwRatio*picturewidth]);
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);

% figure(fig);
% ax = gca;
% exportgraphics(ax, imgName,'Resolution',500);

export_fig(imgName, fig, '-r500', '-transparent');

% produces images of different size when plot data changes but axes remain
% the same
% exportgraphics(ax, imgName,'Resolution',300);
end