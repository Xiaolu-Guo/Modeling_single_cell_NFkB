% 3cm : 85pts
% 3.5cm 100pts
% 4cm : 113pts
% 4.5cm : 128pts

% fontsize: 6-8pts
% title 9-10pts
% abstract 12pts

% medium size figure
paperpos = [0,0,120,100];
papersize = [120,100];
draw_pos = [20,20,70,40];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
