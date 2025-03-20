% Your data matrix
% MI for TAK1, IKK, NFkB, dynamic feature and traj

dataMatrix = [0.0674666787863414, 0.0194968687569173, 0.0305908393000465, 0.0197081010277458, 0.142616265557449;
    0.0505981311504982, 0.0820274472399021, 0.0815022321482828, 0.0680055352468789, 0.0676813601206566;
    0.105062167945661, 0.105062167945718, 0.105062167945718, 0.000235584722696387, 0.121923127044787;
    0.100677780112731, -4.48605703695648e-06, 0.00987832288501522, -1.15920603320774e-05, -8.37114572505016e-05;
    1.04444487256841, 1.05838527679039, 1.06062449088561, 1.05597387788482, 1.05826123609752;
    -0.000621651333801765, -0.000285887582094802, -0.000459500178692451, -0.000374228018245049, -0.000374228018188205];

% Calculating means and standard deviations
means = mean(dataMatrix, 2);
stds = std(dataMatrix, 0, 2);

% Number of categories
numCategories = size(dataMatrix, 1);

X = categorical({'TAK1 dynamic', 'TAK1 time traj', 'IKK dynamic', 'IKK time traj', 'NFkB dynamic', 'NFkB time traj'});
X = reordercats(X,{'TAK1 dynamic', 'TAK1 time traj', 'IKK dynamic', 'IKK time traj', 'NFkB dynamic', 'NFkB time traj'});

% Creating the bar plot
figure;
bar(X,means);

hold on;

% Adding error bars with line width 2
for i = 1:numCategories
    errorbar(i, means(i), stds(i), 'k', 'LineWidth', 2);
end

% Setting labels and title
ylabel('CC');
% title('Bar Plot with Mean and Standard Deviation');

% Ending the hold
hold off;
