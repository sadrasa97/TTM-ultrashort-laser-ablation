%kb=8.6173*10^-5;%ev
kb=1.38*10^-23;
natom=1159;%Number of atoms.
a=3.52*10^-10;
boxsize=35.3*10^-10;%The side-length of the cubic box, unit angstrom.
mass=9.7*10^-26;%Mass of one atom.
%0.007 0.7
%0.2 2.2
temp=300;%Temperature, unit K.
step=100;%Total running step.
epsilon=0.2*1.6*10^-19;%Parameter for LJ potential, the depth of the potential well.(ev)
sigma=2.2*10^-10;%Parameter for LJ potential, the distance at which the inter-particle potential is 0.
rc=2.5*sigma;%Cutoff distance, the LJ potential is truncated at the cutoff distance.
%ecut=4*epsilon*((sigma/rc)^12 - (sigma/rc)^6);%LJ potential at the cutoff distance.
V=35.3*35.3*10*10^-30;%A^3
dt_md=10*10^-15;%fs
%Initilization position.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 142*10*10*4;
%This array will gives is a matrix of Nx3, that will store data.
FCClattice = zeros(N,3);
%Create the vectors for the atoms position in the FCC Lattice:
FCCatoms = [0 0 0;(a/2) (a/2) 0;(a/2) 0 (a/2);0 (a/2) (a/2)];
%Set a variable to change rows in the array storing the coordinates:
 n = 0;
%Create 3 For loops to mix the positions between xyz
for x = 0:10-1
    for y = 0:10-1
        for z = 0:142-1
            for i = 1:4
                %Create a vector that translate your location in the
                %coordinate system into the position in the crystal:
                coordinatestranslation = a*[x y z];
                %Add 1 to move one row down in the array:
                n = n+1;
                FCClattice(n,:) = coordinatestranslation +FCCatoms(i,:);
            end
        end
    end
end

FCClattice=FCClattice';
p=FCClattice(:,1:natom);

%p_ref=p; %Choose initial position as reference for diffusion coefficient calculation.

%Calculate potential and force for the initial condition.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pen=0;%pen is the potential energy.

f=zeros(3,natom);%f is the force.
for i=1:(natom-1)
    for j=(i+1):natom
        r1=p(:,i)-p(:,j);
            r1(1,:)=r1(1,:)-boxsize*round(r1(1,:)./boxsize);%PBC x 
            r1(2,:)=r1(2,:)-boxsize*round(r1(2,:)./boxsize);%PBC y %Consider PBC to get the correct distance between two atoms.
        
        r=sqrt(sum(r1.^2));
        if ((r) <= (rc))
            ff=24*epsilon*sigma^6*(2*sigma^6-r^6)/(r^14);
            f(:,i)=f(:,i)+ff*r1;
            f(:,j)=f(:,j)-ff*r1;
           
            
        end
    end
end
acce=f/mass;%Calculate acceleration on each atom.
aaa=zeros(step,1);
pressure=zeros(step,1);
virial=zeros(step,1);
%Initilization velocity.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=randn(3,natom); %Generate random numbers as initial velocities.
ave_v=sum(v,2)/natom; %Calculate the average value of velocities.
v=(v-repmat(ave_v,[1,natom])); %Scale velocities to make sure the total momentum is zero.
ave_v=sum(v,2)/natom;
ave_v2=sum(sum(v.^2,2)/natom);
fa=sqrt(3*kb*temp/(mass*ave_v2)); %Velocity scale factor to specific temperature.
v=(v-repmat(ave_v,[1,natom]))*fa; %Scale velocites to make sure the total momentum is zero and fit the temperature.
ken=sum(sum(0.5*mass*v.^2,2)); %ken is the kinetic energy.

for m=1:step
    
    p=p+v*dt_md+0.5*(acce+v*kesiprime(m,1))*dt_md^2; %Update position.
 %   p=p+dt_md*(1-(dt_md/2)*(kesiprime(m,1)))*v+f*(dt_md^2/(2*mass));
    %Calculate force and potential according to the new position.
    %pen=0;
   
   
    f=zeros(3,natom);
    for i=1:(natom-1)
        for j=(i+1):natom
            r1=p(:,i)-p(:,j);
            r1(1,:)=r1(1,:)-boxsize*round(r1(1,:)./boxsize);%PBC x 
            r1(2,:)=r1(2,:)-boxsize*round(r1(2,:)./boxsize);%PBC y
            r=sqrt(sum(r1.^2));
            if ((r) <= (rc))
                ff=24*epsilon*sigma^6*(2*sigma^6-r^6)/(r^14);
                f(:,i)=f(:,i)+ff*r1;
                f(:,j)=f(:,j)-ff*r1;

                virial(m,1)=virial(m,1)-ff*r^2;
               
    %            pen=pen+4*epsilon*((sigma./r).^12-(sigma./r).^6)-ecut;
            end
        end
    end


  %  kesiprime=((natom)^-1)*(g*V)*(((Te(m,num_layer)-Ti_md(m,num_layer))+(((Te(2*m,num_layer)))-Ti_md(m,num_layer))))*((3*natom*kb*Ti(m,num_layer))^-1);

    v=(v+0.5*(dt_md*acce+(((f/mass)+((kesiprime(m,1))*v))*dt_md)));%Update velocity.
   
 %  v=(((1-((dt_md/2)*kesiprime(m,1)))/(1+((dt_md/2)*kesiprime(m,1))))*v+(1/(1+((dt_md/2)*kesiprime(m,1))))*(dt_md)*(((acce*mass)+f)/(2*mass)));
 % v=v*((1+((dt_md/1000)*(((Ti(2*m,1))/(mass*sum(sum(v.^2,2)/natom)/(3*kb)))-1)))^0.5);
   
 %  v=v*(sqrt(Ti(8*m,1)/(mass*sum(sum(v.^2,2)/natom)/(3*kb))));
 %  end
    aaa(m,1)=(mass*sum(sum(v.^2,2)/natom)/(3*kb));%Calculate real temperature.
    
    pressure(m,1)=((((natom*kb*(aaa(m,1))/V)))+((virial(m,1))/(3*V)))*10^-9;
   
    
end