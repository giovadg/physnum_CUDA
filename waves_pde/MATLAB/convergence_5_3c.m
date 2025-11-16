% convergence 5.3(c)
nx=[16      32       64        128        256        512 1024];
err=[1.5153 0.655152 0.3166727 0.15682077 0.07818907 0.0390710506 0.01953239];
figure
loglog(nx,err,'k+-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('n_x')
ylabel('err [a.u]')
grid on
%---
% Trick#1 : enlarge domain artificially by "Delta x / 2"
% and choose nx multiple of 2*n+1
% Here n=3
nx=[35 70 140 280 560 1120];
nsteps=[20 40 80 160 320 640];
err=[8.006e-4 4.7875e-5 2.9348e-6 1.81779e-7 1.1312e-8 7.055e-10];
hold on
loglog(nx,err,'r+-','linewidth',lw)

