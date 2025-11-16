% Ex5_2024 5.4d var xa
xb=700e3; hc=35; hl=7000;
xa=[300 600 650 680]*1e3;
Amax=[3.591 3.119 2.675 2.075];
Amax=[(hl/hc)^0.25 Amax];
Arefl=[0.0808 0.3145 0.5746 0.8397];
fs=16; lw=1.5; ms=8;
figure
xinv=[0 1./(xb-xa)];
plot(xinv,Amax,'k+-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('1/(x_b-x_a) [m^{-1}]')
ylabel('Amplitude [m]')
hold on
plot(0,(hl/hc)^0.25,'ko','markersize',ms)
