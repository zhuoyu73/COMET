function [kx, ky, kz]=kspace_add_zencodes(kx_in,ky_in,kz_sample,Kz0Pts)

npts = length(kx_in);
b_xswitch = 0;
b_yswitch = 0;
% kz_step = [0 0.005 0.01 0.02 0.05 0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 0.74 0.81 0.85 0.90 0.93 0.95 0.98 0.99 0.995 1];
tmp = 0:1:15;
tmp = tmp.^2;
tmp = tmp/max(tmp)/2;
tmp_reverse = 1- flipdim(tmp,2);
kz_step = [tmp tmp_reverse(2:end)];

kx = kx_in;
ky = ky_in;
kz = zeros(size(kx_in));
kz_index = 1;
if nargin < 4
    Kz0Pts = 200;
end
for xx = Kz0Pts:npts %don't start doing z steps right away
    if ((sign(kx_in(xx))~=sign(kx_in(xx-1))))&&(kx_in(xx)>kx_in(xx-1))
        b_xswitch = 1;
    end
    if ((sign(ky_in(xx))~=sign(ky_in(xx-1))))&&(ky_in(xx)>ky_in(xx-1))
        b_yswitch = 1;
    end
    if b_xswitch&&b_yswitch
        kz_index_next = mod(kz_index+1,length(kz_sample));
        if kz_index_next == 0
            kz_index_next = length(kz_sample);
        end
        kz((xx-length(kz_step)+1):xx) = kz_sample(kz_index)+kz_step*(kz_sample(kz_index_next)-kz_sample(kz_index));
        b_xswitch = 0;
        b_yswitch = 0;
        kz_index = kz_index_next;
    else
        kz(xx) = kz_sample(kz_index); 
    end
    
end
