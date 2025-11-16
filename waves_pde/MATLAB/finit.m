function [finit]=finit(x,x0,xr,xl,sigma,n,x1,x2,nfunc)
%FINIT fonction f en t=0
[ii,znx]=size(x);
if nfunc==0
   finit=zeros(1,znx); % nulle partout
end
if nfunc==1
   finit=exp(-(x-x0).^2/(2*sigma^2)); % 1:Gaussienne
end
if nfunc==2
   finit=(x-xl).*(xr-x).*exp(-(x-x0).^2/(2*sigma^2)) ...
        / ((x0-xl)*(xr-x0)); % 2: "Gauss" bords nuls
end
if nfunc==3
   finit=sin(2*pi*n*(x-xl)/(xr-xl)); % 3: mode propre bords fixes
end
if nfunc==4
   finit=cos(2*pi*n*(x-xl)/(xr-xl)); % 4: mode propre bords libres
end
if nfunc==5
   zx=xl+mod(x-xl,xr-xl);
   finit=min((zx-xl)/(x0-xl),(xr-zx)/(xr-x0));
%%%   finit=min((x-xl)/(x0-xl),(xr-x)/(xr-x0)); % 5: pincement corde (guitare)
end
if nfunc==6
    finit=zeros(1,znx);
    for j=1:znx
    if ((x(j)>x1)&(x(j)<x2))
        finit(j)=0.5*(1-cos(2*pi*(x(j)-x1)/(x2-x1)));
    end
    end 
end
if nfunc==7 % Mode propre bord fixe a gauche, libre a droite
    finit=sin(pi*(n+0.5)*(x-xl)/(xr-xl)); % 3: mode propre bords fixes
end
    


