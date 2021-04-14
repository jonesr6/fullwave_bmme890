%addpath /mnt/intra_raid6/rebecca/path/
%addpath /mnt/intra_raid6/rebecca/Monkey_skull_CT_sim/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebecca Jones
% WRITTEN: April 12, 2021
% LAST MODIFIED: April 13, 2021
% Simple transcranial fullwave launcher example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;              % speed of sound (m/s)
omega0 = 2*pi*0.5e6;    % center radian frequency of transmitted wave
wY = 6e-2;              % width of simulation field (m)
wZ = 10e-2;             % depth of simulation field (m)
duration = 2.0*wZ/c0;   % duration of simulation (s)
p0 = 1e5;               % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 15;               % number of points per spatial wavelength
cfl = 0.2;              % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
dY = c0/omega0*2*pi/ppw
dZ = c0/omega0*2*pi/ppw
dT = dY/c0*cfl; 
nY=round(wY/dY);            % width of simulation field in pixels
nZ=round(wZ/dZ);            % depth of simulation field in pixels
duration = 2.0*nZ*dZ/c0;    % duration of simulation (s)
nT = round(duration*c0/lambda*ppw/cfl);
%% map conversion %%

%%%%%%%load in skull maps created from CT scans
%load('maps')

%%%%%%%create simple skull model

%%%%%flat skull
%cmap=ones(nY,nZ)*1540;
%skull_start=round(40e-3/dZ);
%skull_end=round(46.5e-3/dZ);
%cmap(:,skull_start:skull_end)=2900;

%%%%%curved skull
cmap=ones(nY,nZ)*1540; %speed of sound in water
ROC=round(.075/dZ); %Radius of curvature of skull model [m]
cen=[round(nY/2) round(115e-3/dZ)]; %center of skull
skull_thickness=round(6.5e-3/dZ);
for ii=1:nY
    for jj=1:nZ
         rr=(ii-cen(1))^2+(jj-cen(2))^2; %calculate radius at each point
         if(rr<(ROC+skull_thickness)^2 && rr>(ROC-0.01)^2)
            cmap(ii,jj)=2900; %speed of sound in dense cortical bone [m/s]
         end
    end
end

figure;imagesc(cmap');

%%%%%%% create remaining maps
rhomap=ones(nY,nZ)*1000; %density of water
rhomap(cmap==2900)=2200; %density of dense cortical bone
imagesc(rhomap')

Amap=zeros(nY,nZ);
phi = 1-mat2gray(cmap); %porosity
idi = find(mat2gray(cmap)>1e-2);
Amap(idi) = (2+(80-2)*sqrt(phi(idi)))*12/40; %attenuation calculation based on porosity
imagesc(Amap');

boveramap = zeros(nY,nZ); % nonlinearity map 

%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nY,nZ); 
inmap(:,1) = ones(nY,1); inmap(:,2) = ones(nY,1); inmap(:,3) = ones(nY,1);
incoords = mapToCoords(inmap);

%%% Generate initial conditions based on input coordinates %%%%%%

% generate focused transmit
foc=round(8e-2/dZ); % location of focus in axial direction
fcen=[round(nY/2) foc]; % center of focus [lateral, axial]

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold all
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold off
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);
imagesc(icmat)

%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap=zeros(nY,nZ); modY=2; modZ=2; 
[modidy modidz]=meshgrid(1:modY:nY,1:modZ:nZ);
outmap(modidy,modidz)=1;
imagesc(outmap');
outcoords=mapToCoords(outmap);

%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd=pwd; addpath(cwd);
outdir=['transcranial/']; eval(['!mkdir -p ' outdir]); 
eval(['!cp try6_nomex ' outdir]);
cd(outdir)
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords,outcoords,icmat);
eval('!./try6_nomex')
cd(cwd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% process sim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
ncoordsout=size(outcoords,1)
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout; nRun=floor(nRun)-1
    
%find pressure
genout = readGenoutSlice([outdir 'genout.dat'],1:nRun,size(outcoords,1));
idc=1:ncoordsout;
p1 = reshape(genout(:,idc),size(genout,1),length(0:modY:nY-1),length(0:modZ:nZ-1));

%pressure movie
for i=1:10:size(p1,1)
     tmp=interp2easy(squeeze(p1(i,:,:)),modZ,modY);
     figure(1);imagesc((1:nY)*dY*100,(1:nZ)*dZ*100,tmp')
     axis equal tight
     i
end

%calculate intensity
pp=squeeze(sum(p1.^2));

%intensity plot
imagesc((1:nY)*dY,(1:nZ)*dZ,squeeze(sum(p1.^2))')

%intensity movie
    tmp=interp2easy(squeeze(p1(1,:,:)),modZ,modY); %interpolates if modY,modZ>1
    pI2x=zeros(size(tmp));
    for i=1:size(px,1);
    pI2x=pI2x+interp2easy(squeeze(p1(i,:,:)),modZ,modY).^2;
    if(~mod(i,10))
        dbz=dbzero(pI2x);
        figure(10);imagesc((1:nY)*dY*100,(1:nZ)*dZ*100,dbz',[-40 0]);colorbar;
        axis equal tight
        i
    end
    end



