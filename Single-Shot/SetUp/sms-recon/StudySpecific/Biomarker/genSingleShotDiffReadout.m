clear
gamp = 22;
gslew = 135;
N = 150;
nl = 2.0;
tsamp = 10e-6;
D = 240;
nsl = 9;
RUP = 40;
RDP = 4;

[Gx, Gy, kx, ky, sx, sy]=genspi(D/10, N, nl, gamp/10, gslew, tsamp);

kz_sample = [0 -4 4];

[kx ky kz] = kspace_add_zencodes(kx(:),ky(:),kz_sample);
% kz = zeros(size(kx));


%[gx gy gz] = kspace2gradients(flipdim(kx(:),1),flipdim(ky(:),1),flipdim(kz(:),1),N,N,nsl, tsamp, RUP, RDP);
grads_mTm = calcGradsFromKSpace([kx(:),ky(:),kz(:)], tsamp, D); 

shotSchedule = [0,0]; %Specific to Single Shot Case

writeExternalGradFileMultiShot('BioMGrads.xml', grads_mTm, shotSchedule,D,1, 0);