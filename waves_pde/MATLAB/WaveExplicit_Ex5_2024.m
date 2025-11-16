% R\'esout l'Eq. d'Alembert (x,t) avec sch\'ema explicite \`a 3 niveaux
% Inclut l'option vitesse de phase variable, u^2=u^2(x)
% Soit u^2(x) d^2f/dx^2 (nopt=1); d/dx(u^2(x)df/dx) (nopt=2); d^2/dx^2(u^2f) (nopt=3)

% nvphase=3; % 0: vitesse de phase uniforme; 
           % 1: coupure en x=xcut, r\'esonance en x=xres; 
           % 2: reef&beach (linear)
           % 3: reef (linear) & beach (power law smooth join)
	   % 4: flat ocean bottom(h0) in [xl,xocean], sin^2 in [xocean,xr], hbeach en x=xr
	   % 5: flat ocean bottom(hocean) in [xl,xrecifl], sin^2 in [xrecifl,xplagel], flat hrecif in [xplagel,xplager], 
	   %    sin^2 in [xplager,xrecifr], flat ocean bottom (hocean) in [xrecifr,xr]
       % 6: Ex.7 2017; 
       % 7: Ex.7 2018; [xr xl]=[-L +L]; h=h1 (x<-r), h2 (x>+r), ~sin(-r<x<+r)
       % 8: Ex.7 2020 (Belharra)
       % 9: Ex.7 2023 ~cos
       % 10: Ex.5 2024 Coral Reef

nopt=2; % 1: u^2(x) d^2f/dx^2; 2: d/dx(u^2(x)df/dx); 3: d^2/dx^2(u^2f)
g=9.81; 
%nx=64; % nx=1024 or 2048:TSUNAMI %%% nx= 64 or 128; % % autres
%------------------------------------------------------------------
% Ex.5 2024 
x0=1; sigma=1; n=1; % dummies
% 5.3 h0=const
nvphase=0;xl=0; xr=10; h0=3; nx=64; A=1; tfin=4*(xr-xl)/sqrt(g*h0);
nbcleft=0; nbcright=1; nfunc=6; x1=2; x2=6; npropag=0; % for 5.3(a)(b)
% tfin=3; % for 5.2(b)
nfunc=7; n=3; tfin=2*(xr-xl)/(sqrt(g*h0)*(n+0.5)); nx=70;% for part 5.3(c)
% Trick#2: calculate tfin as if L--> L-dx/2
% tfin=2*(xr-(xr-xl)/(2*nx))/(sqrt(g*h0)*(n+0.5));
% 5.4 h0 coral reef
% xl=0; xr=1e6; xa=697e3; xb=700e3; xc=720e3; xd=850e3;
% nvphase=10; nx=4096; hl=7000.; hc=35.; hr=200; h0=hl;
% nbcleft=5; nbcright=5; nfunc=6; A=1; x1=50e3; x2=250e3; 
% tfin=12000; npropag=1; T=(x2-x1)/sqrt(g*hl);
% 
% 
%------------------------------------------------------------------
% Ex.7 2023 (~cos)
% nvphase=9; nx=2048;
% hl=7500; hr=20;
% xl=0; xr=1000e3; xa=900e3; xb=950e3; 
% A=1; T=15*60; tfin=10000; %nbcleft=3; nbcright=5; 
% h0=hl;
%------------------------------------------------------------------
% Ex.7 2020 (Belharra)
% nvphase=8; xl=-5000.0; xr=2500.0; hl=750.0; hr=28.0; hhf=230.0; xhf=0.0; sigmahf=500.0;
% A=2.0; T=15.0; tfin=150.0; %nbcleft=3; nbcright=5; 
% h0=hl;
%------------------------------------------------------------------
% Ex.7 2019 (tsunami on Great Barrier reef)
%  nvphase=5; L=800.0e3; xl=0.0; xr=L; h1=8000.0; h2=20.0;  
%  xa=200.0e3; xb=370.0e3; xc=430.0e3; xd=600.0e3;
%  A=1.0; T=900.0; tfin=10000.0; %nbcleft=3; nbcright=5;
%  h0=h1; 
%--------------------------------------------------------------------
% Ex.7 2018 (tsunami, [xr xl]=[-L +L]; h=h1 (x<-r), h2 (x>+r), ~sin(-r<x<+r)
%  nvphase=7; L=500.0e3; xl=-L; xr=L; r=480.0e3; h1=9000.0; h2=20.0; 
%  A=1.0; T=900.0; tfin=12000.0; %nbcleft=3; nbcright=5;
%  h0=h1;
%-------------------------------------------------------------------
% Linear + smooth join to beach
%  nvphase=3; xl=-1000.0e3; xr=0.; x0=xl; h0=7000.; hreef=400.; hbeach=5.;
%  A=1.0; T=900.; tfin=10000; 
%  hreef=400.0; hbeach=20.0; xocean=xl+980.0e3;
% T=900.0; %Ex.7 2017 % T=2000; TSUNAMI % T = 2.0*(xr-xl)/(n*u0); % excitation period
% A=2.0; % TSUNAMI % A=0.1; % MODES PROPRES % % amplitude d'excitation%xl=-1.0e6; xr=0.0; % TSUNAMI % xl=0.0; xr=1.0; % 
%%%tfin=((xr-xl)/u0)*3; 
%tfin=13000.0; % TSUNAMI % tfin=3.53; % TESTS % tfin=5.0; % avec CFL=1.0001, nx=128, INSTABILITY % tfin=10; % MODES PROPRES % 
dx=(xr-xl)/nx;
x=[xl:dx:xr];
% Trick#1: enlarge the domain by "half a grid point"
xrtilde=(xr-xl/(2*nx))/(1-1/(2*nx));
% xrtilde=xr+(xr-xl)/(2*nx); % less accurate
dx=(xrtilde-xl)/nx;
x=[xl:dx:xrtilde];

[ii,NX]=size(x);

plot_time_evolution=0; % 0/1: do not / do plot f(x,.) at each timestep
plot_surf=1; % 0/1 do not / do contourf plot fall(x,t)

% nbcleft=4;  % bord gauche 0:fixe; 1:libre; 2:p\'eriodique; 3: A*sin(om*t) 4: A*sin(om*t) nper periods; 
% 	    %             5: sortie 6: A*sin(om*t)*(1/(3*sqrt(3)/4))*(1-cos(om*t)) nper periods;
% nbcright=5; % bord droite 0:fixe; 1:libre; 2:p\'eriodique; 5: sortie
% 
% npropag=1;  % propagation -1: vers la gauche; 0: gauche et droite; +1: droite
% 		   %  si u positif; autrement gauche et droite invers\'es!
% 
% nfunc=0; % fonction initiale 0: nulle; 1: Gauss; 2: "Gauss" bords nuls; 3: sin; 4: cos; 5: pincement corde

CFL=1.000; % >1.0: instable!!!  u^2(x)<0 instable!!!

% Ex.7 2014
%hocean=6000.; hrecif=25.;
%%%xl=-1000.0e3; xr=0.;<
%xrecifl=-700.0e3; xplagel=-530.0e3; xplager=-470.0e3; xrecifr=-300.0e3;
%A=1.0; 
%h0=hocean;

u0=sqrt(g*h0); % TSUNAMI % u0=1.0;  %  vitesse de phase % 

if (nvphase==2 | nvphase==3)
    ureef=sqrt(g*hreef); % vphase pr\`es du r\'ecif (nvphase=2)
    ubeach=sqrt(g*hbeach); %vphase pr\`es de la plage
    xreef=xl+(xr-xl)*0.8; % position du "r\'ecif"
    xbeach=xr+(xr-xl)*0.00; % position de la "plage"
end

dt=CFL*dx/abs(u0)  
if nfunc==7
    nsteps=ceil(4*nx/(2*n+1)); dt=tfin/nsteps; CFL=u0*dt/dx;
end
t=[0:dt:tfin];
[ii,nt]=size(t)

% x0 = xl + (xr-xl)*0.1; % centre Gaussienne initiale ou pincement initial
% sigma = (xr-xl)*0.025;  % ecart-type de la Gaussienne initiale
% n=4.0; %  nb de demi ong.d'ondes (=multiple de la fondamentale  du mode propre \`a birds fixes)
%        % int\'eressant: entier ou demi-entier, nbcright=0 or =1  
% om=2*pi/T; %:TSUNAMI %%% om=n*pi*u0/(xr-xl); %: EIGENMODES
% nper=1; % dur\'ee de l'excitation en nombre de p\'eriodes (option nbcleft=4 ou 6)


% xres=xl+(xr-xl)*1.001;
% xcut=xl+(xr-xl)*1.5;
r2=zeros(1,NX);
u2=zeros(1,NX);
if nvphase==0
   u2=u0^2*ones(1,NX);
elseif nvphase==1
   u2=u0^2*(x-xres)./(x-xcut);
elseif nvphase==2
   for jx=1:NX
      if (x(jx)<xocean)
         u2(jx)=u0^2+(x(jx)-xl)*(ureef^2-u0^2)/(xreef-xl);
      else
         u2(jx)=ureef^2+(x(jx)-xreef)*(ubeach^2-ureef^2)/(xbeach-xreef);
      end
   end
elseif nvphase==3
   zs=(ureef^2-u0^2)/(xreef-xl); % slope in deep ocean
   zbeta=-zs*(xbeach-xreef)/(ureef^2-ubeach^2);
   zalpha=-zs/(zbeta*(xbeach-xreef)^(zbeta-1));
   for jx=1:NX
      if (x(jx)<xreef)
         u2(jx)=u0^2+(x(jx)-xl)*(ureef^2-u0^2)/(xreef-xl);
      else
         u2(jx)= ubeach^2 + zalpha*(xbeach-x(jx))^zbeta;
      end
   end
elseif nvphase==4
   for jx=1:NX
      if (x(jx) < xocean)
	 u2(jx) = g*h0;
      else
	 u2(jx) = g*h0 + g*(hbeach-h0)*(sin(pi*(x(jx)-xocean)/(2.0*(xr-xocean))))^2;
      end
   end
elseif nvphase==5
   for jx=1:NX
      if (x(jx) < xa)
	 u2(jx) = g*h1;
      elseif (x(jx) < xb)
	 u2(jx) = g*h1 + g*(h2-h1)*(sin(pi*(x(jx)-xa)/(2.0*(xb-xa))))^2;
      elseif (x(jx) < xc)
	 u2(jx) = g*h2;
      elseif (x(jx) < xd)
	 u2(jx) = g*h1 + g*(h2-h1)*(sin(pi*(x(jx)-xd)/(2.0*(xc-xd))))^2;
      else
	 u2(jx) = g*h1;
      end
   end
elseif nvphase==6
   for jx=1:NX
      if (x(jx) < xocean)
         u2(jx) = g*(h0+(hbeach-h0)*(sin(0.5*pi*(x(jx)-xl)/(xocean-xl))^2));
      else
         u2(jx) = g*hbeach;
      end
   end
elseif nvphase==7
   for jx=1:NX
      if (x(jx) < -r)
        u2(jx) = g*h1;
      elseif (x(jx) < r)
%         u2(jx) = g*(h1+(h2-h1)*(sin(0.5*pi*(x(jx)+r)/(2.0*r))^2));
        u2(jx) = g*0.5*((h1+h2)+(h2-h1)*(sin(0.5*pi*x(jx)/r)));
      else
        u2(jx) = g*h2;
      end
   end
elseif nvphase==8
    for jx=1:NX
        u2(jx)=g*(hl+(hr-hl)*(x(jx)-xl)/(xr-xl) - hhf*exp(-(x(jx)-xhf)^4/sigmahf^4));
    end
elseif nvphase==9
    for jx=1:NX
        if x(jx)<=xa
            u2(jx)=g*hl;
        elseif x(jx)<=xb
            u2(jx)=0.5*g*(hl+hr + (hl-hr)*cos(pi*(x(jx)-xa)/(xb-xa)));
        else
            u2(jx)=g*hr;
        end
    end
elseif nvphase==10
    for jx=1:NX
        if x(jx)<=xa
            u2(jx)=g*hl;
        elseif x(jx)<xb
            u2(jx)=0.5*g*(hl+hc + (hl-hc)*cos(pi*(x(jx)-xa)/(xb-xa)));
        elseif x(jx)<=xc
            u2(jx)=g*hc;
        elseif x(jx)<xd
            u2(jx)=0.5*g*(hr+hc - (hr-hc)*cos(pi*(x(jx)-xc)/(xd-xc)));
        else
            u2(jx)=g*hr;
        end
    end
end

r2=u2*dt^2/dx^2;

fm=zeros(1,NX); % f au pas temporel pr\'ec\'edent
f=zeros(1,NX); % f au pas temporel actuel
fp=zeros(1,NX); % f au pas temporel suivant

fall=zeros(nt,NX); % f(t,x)

f=finit(x,x0,xr,xl,sigma,n,x1,x2,nfunc);

fall(1,:)=f;

if npropag==0
   fm=f; % t at timestep -1 => propagation gauche et droite
else
% propagation initiale vers la droite ou vers la gauche (seulement cas u uniforme!!!)
% vers la droite
   if npropag==1
      xsh=u0*dt;
   end
% vers la gauche
   if npropag==-1
      xsh=-u0*dt;
   end
   fm=finit(x+xsh,x0,xr,xl,sigma,n,x1,x2,nfunc);
end

fmax=max(f); fmin=min(f); 
if nvphase==10
    fampli_init=fmax; % positive height of the wave
else
    fampli_init=0.5*(fmax-fmin);
end
if nfunc==0
   fampli_init=abs(A);
end
fabsmax=(max(u2)/min(u2))^0.25*abs(A);

str=('Onde Explicite 3 niveaux')
fs=16; lw=1.5;

%if (nvphase~=1)
   figure
   h=plot(x,-u2/g,'k-','linewidth',lw);
   set(gca,'fontsize',fs)
   xlabel('x [m]')
   ylabel('z [m]')
   title([str])
%end

figure
hf=plot(x,f,'bo-', x,fm,'r*-');
set(gca,'fontsize',fs)
xlabel('x [m]')
ylabel('f [a.u.]')
title([str])
zfactor=1.2; % 
axis([xl xr -zfactor*fabsmax +zfactor*fabsmax]) 
%axis('manual')
%%%%%%M(1)=getframe %%% careful with movies: 32 intervals, 129 timesteps ok
pause
xtext=xl+0.05*(xr-xl);
ytext=0.1*fabsmax;

% timesteps ------------------------------------------------------------------------------

for jt=1:nt-1

% boucle sur x ---------------------------------------------------------------------------
   for jx=2:NX-1
      if (nopt==1 | nopt==2) 
	fp(jx) = 2*(1-r2(jx))*f(jx) - fm(jx) + r2(jx)*(f(jx+1)+f(jx-1));
	  if nopt==2
	   fp(jx) = fp(jx) + 0.25*(r2(jx+1)-r2(jx-1))*(f(jx+1)-f(jx-1));
	  end
      elseif (nopt==3)
	fp(jx) = 2*(1-r2(jx))*f(jx) - fm(jx) + r2(jx+1)*f(jx+1) + r2(jx-1)*f(jx-1);
      end
   end
   
% avec des conditions aux bords fixes (Dirichlet), 
% ne pas changer f en jx=1 ni en jx=NX
   if nbcleft==0
      fp(1)=f(1);
   end
   if nbcright==0
      fp(NX)=f(NX);
   end

% avec des conditions bords libres (Neumann),
% la d\'eriv\'ee est nulle
   if nbcleft==1 
      fp(1)=fp(2);
   end
   if nbcright==1
      fp(NX)=fp(NX-1);
   end

% avec des conditions aux bords p\'eriodiques
   if nbcleft==2
      fp(1) =  2*(1-r2(1))*f(1) - fm(1) + r2(1)*(f(2)+f(NX-1));
   end
   if nbcright==2
      fp(NX) =  2*(1-r2(NX))*f(NX) - fm(NX) + r2(NX)*(f(2)+f(NX-1));
   end

% avec excitation au bord gauche
   if nbcleft==3
      fp(1)=A*sin(om*dt*jt);
   end
 
% avec excitation au bord gauche pendant nper p\'eriodes
   if nbcleft==4
      if (dt*(jt-0.5)<nper*2*pi/om)
         fp(1)=A*sin(om*dt*jt); 
      elseif (dt*(jt-1.5)<nper*2*pi/om)
         fp(1)= 0.0;
      else
        fp(1) = f(1) + (f(2)-f(1))*sqrt(r2(1)); %%%%0.0;
      end
   end
   
% onde sortant par la gauche
   if nbcleft==5
      fp(1) = f(1) + (f(2)-f(1))*sqrt(r2(1));
   end

% avec excitation au bord gauche pendant nper p\'eriodes
   if nbcleft==6
      if (dt*jt<nper*2*pi/om)
	fp(1)=sin(om*dt*jt)*(A/(3*sqrt(3)/4))*(1-cos(om*dt*jt)); %%%%%%%%%%%%%%;
      else
         fp(1)=0.0;
      end
   end


% onde sortant par la droite
   if nbcright==5
      fp(NX) = f(NX) -(f(NX)-f(NX-1))*sqrt(r2(NX));
   end
      
%%%   if ((mod(jt,12))==1)&(jt>212)
%%%     hf=plot(x,f-jt/12,'b-');
%%%     hold on
%%%   end
if plot_time_evolution==1
   hf=plot(x,fp,'b-', 'linewidth',lw);
   ht=text(xtext,ytext,num2str(jt*dt,3));
   axis([xl xr -zfactor*fabsmax +zfactor*fabsmax]) 
   xlabel('x [m]')
   ylabel('f [a.u.]')
   title([str,' t=',num2str(jt*dt,3)])
   pause(0.001)
end

%%%%%%   M(jt+1)=getframe;

%stocker f dans le tableau fall
   fall(jt+1,:)=fp;

% pr\'eparer le prochain pas temporel
   fm=f;
   f=fp;

end


%%%%%%movie(M,1)

% plot(s) final(s)
hf=plot(x,fp,'b-', 'linewidth',lw);
ht=text(xtext,ytext,num2str(jt*dt,3));
axis([xl xr -zfactor*fabsmax +zfactor*fabsmax]) 
xlabel('x [m]')
ylabel('f [a.u.]')
title([str,' t=',num2str(jt*dt,3)])

% for nfunc=7 (eigenmode), compute error after 1 period
if nfunc==7
    err=0;
    for jx=1:NX
        err=err+abs(fall(end,jx)-fall(1,jx))*dx;
    end
    err
end
       
%----------------------------------------------------------------
% non-uniform cases
%----------------------------------------------------------------
if not(nvphase==0)
    if nvphase==10
        fampli=max(fall);
    else
        fampli=(max(fall)-min(fall))*0.5; % amplitude as function of x
    end
   vphaserel=sqrt(u2)/u0;
% for initial condition npropag==0, amplitude is halved
   if npropag==0
      fampli_init=fampli_init*0.5;
   end
   figure
   if nopt==1
      hwkb=plot(x,fampli,'k-',x,fampli_init*sqrt(vphaserel),'r-');
   elseif nopt==2
      hwkb=plot(x,fampli,'k-',x,fampli_init./sqrt(vphaserel),'r-');
   elseif nopt==3
      hwkb=plot(x,fampli,'k-',x,fampli_init*vphaserel.^(-1.5),'r-');
   end
   set(hwkb, 'linewidth',lw)
   set(gca,'fontsize',fs)
   xlabel('x [m]')
   ylabel('Amplitude')
%
% Propagation velocity from x positions of peak perturbation at t=t_j
% fall contains f(t_j,x_i)
% Vary the nptsfit and k values until you get "nice" results...
% nptsfit=1, k=20 gives "nice" results
% NB: this technique is only valid for purely propagating waves in a  
% given direction; it becomes wrong as soon as there are reflections!
%
   jstart=1; jend=nt; 
   xcrete=zeros(1,jend-jstart+1);
   tcrete=zeros(1,jend-jstart+1);
   jcount=0; 
   nptsfit=1; % fit using x grid points between ifmax-nptsfit and ifmax+nptsfit 
   [fmax,ifmax]=max(fall,[],2); % maxima of f along x for each time t_j
   for j=jstart:jend
      i=ifmax(j); % index along x of max(f) at fixed t_j
      if ( (i>nptsfit) && (i<NX-nptsfit)) % make the fit only if there are enough grid points available
        jcount=jcount+1;
        xarr=[x(i-nptsfit:i+nptsfit)]-x(i); % data points used for fit
        farr=[fall(j,i-nptsfit:i+nptsfit)]; % data points used for fit
        p=polyfit(xarr,farr,2);             % parabolic fit
        xcrete(jcount)=-p(2)/(2*p(1))+x(i); % position of wave crest
        tcrete(jcount)=t(j);                % time of wave crest
      end
   end

   figure
   plot(tcrete,xcrete,'b+-', 'linewidth',lw)
   set(gca,'fontsize',fs)
   xlabel('t [s]')
   ylabel('x [m]')
   
   k=1; ncrete=length(xcrete); velocity=zeros(1,ncrete);
   for i=1+k:ncrete-k
       velocity(i)=(xcrete(i+k)-xcrete(i-k))/(tcrete(i+k)-tcrete(i-k));
   end

   figure
   plot(xcrete(1+k:ncrete-k),velocity(1+k:ncrete-k),'b+-',x,sqrt(u2),'k--');
   set(gca,'fontsize',fs)
   ylabel('u [m/s]')
   xlabel('x [m]')

end

if plot_surf==1 
    figure
    [X,TGRID]=meshgrid(x,t);
    %zv=[-8:1:8];
    contourf(X,TGRID,fall,15);
    shading flat
    set(gca,'fontsize',fs)
    xlabel('x [m]')
    ylabel('t [s]')
end

% TSUNAMI: comment the following lines
%figure
%hs=surfc(X,TGRID,fall);
%set(gca,'fontsize',fs)
%xlabel('x')
%ylabel('t')
%zlabel('f')
