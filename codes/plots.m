% plot E and M

load('H3AF_E_112');
A = H3AF_E_112;
clear('H3AF_E_112');

plot(A(:,1),A(:,2:size(A,2)));

axis([0 1 -1 -0.8]);
title('N=4096');