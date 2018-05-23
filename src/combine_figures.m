% draws figures on a single graph
% This code could be useful for putting mutliple plots on a single figure
% in the future

path = '../out/synth_cross_corr/small elong/';

% loads figures
a1 = hgload([path 'central_moment.fig']);
a2 = hgload([path 'central_moment_1st_deriv.fig']);
a3 = hgload([path 'central_moment_2nd_deriv.fig']);

%creates figure of subplots
figure;
h(1) = subplot(1,3,1);
h(2) = subplot(1,3,2);
h(3) = subplot(1,3,3);

%loads the invidiuals figures into the combined one

copyobj(allchild(get(a1, 'CurrentAxes')), h(1));
g1 = gca(a1);
title(h(1), g1.Title.String);
copyobj(allchild(get(a2, 'CurrentAxes')), h(2));
g2 = gca(a2);
title(h(2), g2.Title.String);
copyobj(allchild(get(a3, 'CurrentAxes')), h(3));
g3 = gca(a3);
title(h(3), g3.Title.String);

close(a1)
close(a2)
close(a3)
