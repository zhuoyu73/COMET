clear all
gamp = 21;
gslew = 130;
N = 120;
nl_design = 30;
nl_use = 3;
tsamp = 5e-6;
D = 240;
nsl_design = 3;
nsl_use = 3;
RUP = 20;
RDP = 20;
VD = 0;
alphavd = 4;

if VD==1
[Gx, Gy, kx, ky, sx, sy]=genspivd_Kim2(D/10, N, nl_design, gamp/10, gslew, tsamp,alphavd);
else
[Gx, Gy, kx, ky, sx, sy]=genspi(D/10, N, nl_design, gamp/10, gslew, tsamp);
end


phi = 2*pi/nl_use;
for ii = 0:(nl_use-1)
    ang_rot = phi*(ii-(nl_use-1)*floor(ii/nl_use));
    kx_new(:,ii+1) = kx*cos(ang_rot) + ky*sin(ang_rot);
    ky_new(:,ii+1) = ky*cos(ang_rot) - kx*sin(ang_rot);
end

kx_all = repmat(kx_new, [1 nsl_use]);
ky_all = repmat(ky_new, [1 nsl_use]);
kz_all = zeros(size(kx_all));

for ii = 1:nsl_use
    ki1 = (nl_use*(ii-1))+1;
    ki2 = (nl_use*(ii-1))+nl_use;
    
    kz_all(:,ki1:ki2) = ((nsl_design/2)*-1)+(ii-1)*(nsl_design/nsl_use);
end
    
us_mask = 1:size(kx_all,2);
% us_mask = [1 3 10 12];

[gx gy gz] = kspace2gradients(kx_all(:,us_mask),ky_all(:,us_mask),kz_all(:,us_mask),N,N,nsl_design, tsamp, RUP, RDP);


