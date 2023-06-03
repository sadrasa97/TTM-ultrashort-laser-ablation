% TTM-MD code for ultrashort laser ablation of nickel base on D.ianov
clc;
ln2=0.6931471806;
kesiprime=zeros(4000,1);
%F=1.0645*I0*tau
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%laser & target parameter
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
F=430;
%F=430;
I0=5.6493*10^17;
Tm=1439;
T0=300;
tau=200*10^-15;%s
g=3.6*10^17;%electron-phonon coupling coefficient (w/m^3k)
alpha=7.4*10^7;%m^-1
R=0.6;
A=(1-R);%(1-0.6);
ce=1065;%(J/(m^3k^2))
ke0=91;%(w/mk)
kl=97.5;
ci=4.2*10^6;%(J/(m^3K))
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%space cordinate
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

zmax=50*10^-9;%m
dz=1*10^-9;%m
zmin=0;
nz=((zmax)/dz);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%time
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

tmax=50*10^-12;%ps
dt=1.25*10^-15;%ps
tmin=0;
nt=((tmax)/dt)+1;

z=linspace(0,zmax,nz);
z=z';
t=linspace(0,tmax,nt);

%cor=zeros(10000,1000);
%cor(1,:)=z;
%cor(:,1)=t;
%##########################################################################
%initial condition
%rows are time
%##########################################################################
Te=zeros(40000,50);
Te(:,:)=300;

%Te(1,:)=z;
%Te(:,1)=t;
%##########################################################################
Ti=zeros(40000,50);
Ti(:,:)=300;

%Ti(1,:)=z;
%Ti(:,1)=t;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%finite difference
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s=zeros(40000,50);

for n=1:nt
    
   %space 
    for r=1:nz
        %s(n,r)=A*alpha*(I0*exp(-4*ln2*((dt*n)^2)/tau^2)*exp(-alpha*(dz*r)));
         s(n,r)=A*alpha*I0*(exp(-(4*ln2*(t(n))^2)/(tau^2))*exp(-alpha*z(r)));
         
    
         
        
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%boundary condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:(nt-1)
        Te(n+1,1)=Te(n,1)-g*(dt/(ce*Te(n,1)))*(Te(n,1)-Ti(n,1))+(dt/(ce*Te(n,1)))*s(n,1);%z=0 hame zaman ha moshtagh sefr ast.
         
        Te(n+1,50)=Te(n,50)+(dt/(ce*Te(n,50)))*((ke0*(Te(n,50)/Ti(n,50)))^-1)*(1-R)*I0*(exp((-4*ln2*(t(n))^2)/(tau^2)))-g*(dt/(ce*Te(n,50)))*(Te(n,50)-Ti(n,50))+(dt/(ce*Te(n,50)))*s(n,50);

     Ti(n+1,1)=Ti(n,1)+(dt/(ci*3))*(g*(Te(n,1)-Ti(n,1)));

end



%time

%Tf for nickel is 136000k
for n=1:(nt-1)
    
   %space 
for r=2:(nz-1)
     
     Te(n+1,r)=Te(n,r)+(dt/(ce*Te(n,r)))*ke0*(((Te(n,r)/Ti(n,r))*((Te(n,r+1)-2*Te(n,r)+Te(n,r-1))/dz^2)))-g*(dt/(ce*Te(n,r)))*(Te(n,r)-Ti(n,r))+(dt/(ce*Te(n,r)))*s(n,r);

     Ti(n+1,r)=Ti(n,r)+((kl*dt)/ci)*((Ti(n,r+1)-2*Ti(n,r)+Ti(n,r-1))/dz^2)+(dt/ci)*(g*(Te(n,r)-Ti(n,r)));
     
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ti_md=zeros(40000,50);
Ti_md(:,:)=300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nt-1
    Te(n+1,50)=Te(n,50)-(dt/(ce*Te(n,50)))*((ke0*(Te(n,50)/Ti(n,50)))^-1)*(1-R)*I0*(exp((-4*ln2*(t(n))^2)/(tau^2)))-g*(dt/(ce*Te(n,50)))*(Te(n,50)-Ti(n,50))+(dt/(ce*Te(n,50)))*s(n,50);
    
end

natom=1159;
num_layer=1;
V=35.3*35.3*10*10^-30;
kb=1.38*10^-23;

for m=1:5000

%kesiplus(m,1)=((natom)^-1)*(g*V)*(((Te(m,num_layer)-Ti_md(m,num_layer))+(((Te(2*m,num_layer)))-Ti_md(m,num_layer))))*((3*natom*kb*Ti(m,num_layer))^-1);
kesiprime(m,1)=(10/1.25)*(0.28*g*V)*(((Te(m,num_layer)-Ti(m,num_layer))+...
    (((Te(2*m,num_layer)))-Ti(m,num_layer)))+(((Te(3*m,num_layer)))-Ti(m,num_layer))+...
    (((Te(4*m,num_layer)))-Ti(m,num_layer))+(((Te(5*m,num_layer)))-Ti(m,num_layer))+...
    (((Te(6*m,num_layer)))-Ti(m,num_layer))+(((Te(7*m,num_layer)))-Ti(m,num_layer))+...
    (((Te(8*m,num_layer)))-Ti(m,num_layer)))*((3*natom*kb*(sum(Ti(m:8*m,num_layer))/8))^-1);
    %kesitherd(m,1)=((natom)^-1)*(g*V)*(((Te(m,num_layer)-Ti_md(m,num_layer))+(((Te(2*m,num_layer)))-Ti_md(m,num_layer))))*((3*natom*kb*Ti(m,num_layer))^-1);
end
