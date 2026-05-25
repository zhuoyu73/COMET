function [kspace, grads_mTm] = gen_nav_readout(N,multibandFactor, Rxy,Rz, FOV,FOVz, gAmp, gSlew, tSamp)

% gamp = 22;
% gslew = 135;
% N = 150;
% nl = 2.0;
% tsamp = 10e-6;
% D = 240;
gts = 10E-6;
RUP = 50;
RDP = 10;
nShots = Rxy;

[Gx, Gy, kx, ky, sx, sy]=genspi(FOV/10, N, nShots, gAmp/10, gSlew, tSamp);


%Need to design zero connectors for sanity here.
% nPtsRequired = max(ceil(Gx(end)*10/gSlew/gts/1E3), ceil(Gy(end)*10/gSlew/gts/1E3));
% 
% GrampDownX = linspace(Gx(end),0,nPtsRequired);
% GrampDownY = linspace(Gy(end),0,nPtsRequired);
% 
% Gx = [col(Gx); col(GrampDownX)];
% Gy = [col(Gy); col(GrampDownY)];

if multibandFactor > 1
    kz = (-multibandFactor/2):(multibandFactor/2-1);
end
nKz = ceil(length(kz)/Rz);

kzSample = 0;
deltaKz = multibandFactor/nKz;
for ii = 1:(nKz-1)
    PlusMinus = (-1)^ii;
    kzSample(ii+1) = kzSample(ii) + PlusMinus*(deltaKz*ii);
end

[kx, ky, kz] = kspace_add_zencodes(kx(:),ky(:),kzSample,floor(length(kx)/20));
% kz = zeros(size(kx));

% go to GE world;

[k,g,s,m1] = calcgradinfo([Gx(:),Gy(:)],gts,[0,0]);
gRew = mintimegrad(40,[g(end,:),0],[0 0 0],[-k(end,:),-kz(end)/FOV*10],gts,3,gSlew*100,100*gts,3);

%gRew = flip(gRew,1);

%g = vertcat(g,gRew);
[kRew,g,s,m1] = calcgradinfo(gRew,gts,[k(end,:),kz(end)/FOV*10]);

% Get out of GE world 
if(length(k(:,1)) < kz)
    kx = vertcat(0,col(k(:,1))*FOV/10,col(kRew(:,1))*FOV/10);
    ky = vertcat(0,col(k(:,2))*FOV/10,col(kRew(:,2))*FOV/10);
    kz = vertcat(0,zeros(length(col(k(:,1))) - length(kz(:)),1), kz,col(kRew(:,3))*FOV/10);
else
    kx = vertcat(0,zeros(abs(length(col(k(:,1))) - length(kx(:))),1),col(k(:,1))*FOV/10,col(kRew(:,1))*FOV/10);
    ky = vertcat(0,zeros(abs(length(col(k(:,2))) - length(ky(:))),1),col(k(:,2))*FOV/10,col(kRew(:,2))*FOV/10);
    kz = vertcat(0,kz,col(kRew(:,3))*FOV/10);
end

kspace = [flipud(kx), flipud(ky), flipud(kz)];

grads_mTm = calcGradsFromKSpace(kspace,gts,FOV);
grads_mTmZ = calcGradsFromKSpace(kspace,gts,FOVz);

grads_mTm(:,1:2) = grads_mTm(:,1:2);
grads_mTm(:,3) = grads_mTmZ(:,3);


centerOfKSpace = length(kspace(:,1));
filename = ['Rxy' num2str(multibandFactor/nShots) 'Rz' num2str(Rz) 'Nxy' num2str(N) 'Nz' num2str(multibandFactor) '.xml'];
writeExternalGradFile(filename, grads_mTm, FOV, FOVz, centerOfKSpace,N);
save kspace.mat kspace grads_mTm


