function [gx gy gz] = kspace2gradients(kx,ky,kz,Nx,Ny,Nz,tsamp,gAmp, lRampUpPts, lRampDownPts)

%Gref = 26;				% reference gradient in mT/m
FOVref = 240; %in mm
%maxGrad = 26;   %for read phase slice perscription, mT/m
%maxSlew = 210;
ladc = 1024;
%tsamp = 5; %in micorseconds
gts = 10E-6;
gambar = 42.58; 
tCenter = 0;
% FOVminX = 500;  %in mm
% FOVminY = 500;  %in mm
% FOVminZ = 500;  %in mm
lshots = size(kx,2);

%create prepgradient
gx_pre = zeros(lRampUpPts,lshots);
gy_pre = zeros(lRampUpPts,lshots);
gz_pre = zeros(lRampUpPts,lshots);
kx_pre = zeros(1,lshots);
ky_pre = zeros(1,lshots);
kz_pre = zeros(1,lshots);

for ll=1:lshots
    if kx(1,ll) == 0
        gx_pre(:,ll) = 0;
    else
        if (max(kx(:,ll))==min(kx(:,ll))) %readout is only a phase encode
            Gpeak = kx(1,ll)/FOVref/gambar/(lRampUpPts-1)/(gts)*2;
            if mod(lRampUpPts,2)== 0
                Gpeak = Gpeak*(1+1/(lRampUpPts*(lRampUpPts-2)));
            end
            for xx = 1:(lRampUpPts)
                gx_pre(xx,ll) = (1-abs(xx-(lRampUpPts+1)/2)/(((lRampUpPts+1))/2-1))*Gpeak;
            end
        else %we need a ramp to the correct kspace point
            gx_pre(lRampUpPts,ll) = diff(kx(1:2,ll))/gambar/FOVref/tsamp;
            pts_phase_blip = floor(lRampUpPts*2/3);
            pts_for_ramp = floor(lRampUpPts/3);
            for xx = 1:pts_for_ramp
                gx_pre(lRampUpPts+1-xx,ll) = gx_pre(lRampUpPts,ll)/(pts_for_ramp-1)*(pts_for_ramp-xx);
            end
            Garea = kx(1,ll)/FOVref/gambar/(gts)-sum(gx_pre((lRampUpPts-pts_for_ramp):end,ll))-gx_pre(lRampUpPts,ll)*(tsamp/gts);
            Gpeak = Garea/(pts_phase_blip-1)*2;
            if mod(pts_phase_blip,2)== 0
                Gpeak = Gpeak*(1+1/(pts_phase_blip*(pts_phase_blip-2)));
            end
            for xx = 1:pts_phase_blip
                gx_pre(xx,ll) = (1-abs(xx-(pts_phase_blip+1)/2)/((pts_phase_blip+1)/2-1))*Gpeak;
            end
        end
    end
    kx_pre(1,ll) = sum([0; gx_pre(:,ll)])*gts*FOVref*gambar;
    
    if ky(1,ll) == 0
        gy_pre(:,ll) = 0;
    else
        if (max(ky(:,ll))==min(ky(:,ll))) %readout is only a phase encode
            Gpeak = ky(1,ll)/FOVref/gambar/(lRampUpPts-1)/(gts)*2;
            if mod(lRampUpPts,2)== 0
                Gpeak = Gpeak*(1+1/(lRampUpPts*(lRampUpPts-2)));
            end
            for xx = 1:(lRampUpPts)
                gy_pre(xx,ll) = (1-abs(xx-(lRampUpPts+1)/2)/(((lRampUpPts+1))/2-1))*Gpeak;
            end
         else %we need a ramp to the correct kspace point
            gy_pre(lRampUpPts,ll) = diff(ky(1:2,ll))/gambar/FOVref/tsamp;
            pts_phase_blip = floor(lRampUpPts*2/3);
            pts_for_ramp = floor(lRampUpPts/3);
            for xx = 1:pts_for_ramp
                gy_pre(lRampUpPts+1-xx,ll) = gy_pre(lRampUpPts,ll)/(pts_for_ramp-1)*(pts_for_ramp-xx);
            end
            Garea = ky(1,ll)/FOVref/gambar/(gts)-sum(gy_pre((lRampUpPts-pts_for_ramp):end,ll))-gy_pre(lRampUpPts,ll)*(tsamp/gts);
            Gpeak = Garea/(pts_phase_blip-1)*2;
            if mod(pts_phase_blip,2)== 0
                Gpeak = Gpeak*(1+1/(pts_phase_blip*(pts_phase_blip-2)));
            end
            for xx = 1:pts_phase_blip
                gy_pre(xx,ll) = (1-abs(xx-(pts_phase_blip+1)/2)/((pts_phase_blip+1)/2-1))*Gpeak;
            end
        end
    end
    ky_pre(1,ll) = sum([0; gy_pre(:,ll)])*gts*FOVref*gambar;
    
   if kz(1,ll) == 0
        gz_pre(:,ll) = 0;
        kz_pre(:,ll)= 0;
    else
        if (max(kz(:,ll))==min(kz(:,ll))) %readout is only a phase encode
            Gpeak = kz(1,ll)/FOVref/gambar/(lRampUpPts-1)/(gts)*2;
            if mod(lRampUpPts,2)== 0
                Gpeak = Gpeak*(1+1/(lRampUpPts*(lRampUpPts-2)));
            end
            for xx = 1:(lRampUpPts)
                gz_pre(xx,ll) = (1-abs(xx-(lRampUpPts+1)/2)/(((lRampUpPts+1))/2-1))*Gpeak;
            end
            kz_pre(:,ll)= 0;
        else %we need a ramp to the correct kspace point
            gz_pre(lRampUpPts,ll) = diff(kz(1:2,ll))/gambar/FOVref/tsamp;
            pts_phase_blip = floor(lRampUpPts*2/3);
            pts_for_ramp = floor(lRampUpPts/3);
            for xx = 1:pts_for_ramp
                gz_pre(lRampUpPts+1-xx,ll) = gz_pre(lRampUpPts,ll)/(pts_for_ramp-1)*(pts_for_ramp-xx);
            end
            Garea = kz(1,ll)/FOVref/gambar/(gts)-sum(gz_pre((lRampUpPts-pts_for_ramp):end,ll))-gz_pre(lRampUpPts,ll)*(tsamp/gts);
            Gpeak = Garea/(pts_phase_blip-1)*2;
            if mod(pts_phase_blip,2)== 0
                Gpeak = Gpeak*(1+1/(pts_phase_blip*(pts_phase_blip-2)));
            end
            for xx = 1:pts_phase_blip
                gz_pre(xx,ll) = (1-abs(xx-(pts_phase_blip+1)/2)/((pts_phase_blip+1)/2-1))*Gpeak;
            end
        end
   end
    kz_pre(1,ll) = sum([0; gz_pre(:,ll)])*gts*FOVref*gambar;
end

%convert kspace to gradients;
gx = diff([kx_pre(:,:); kx])/gambar/FOVref/tsamp;
gy = diff([ky_pre(:,:); ky])/gambar/FOVref/tsamp;
gz = diff([kz_pre(:,:); kz])/gambar/FOVref/tsamp;
%gy=[gy_pre(lRampUpPts,:); diff(ky)/gambar/FOVref/tsamp*1e6];
%gz=[gz_pre(lRampUpPts,:); diff(kz)/gambar/FOVref/tsamp*1e6];

%interpolate to gradient raster time
tGrad = tsamp*length(kx)+mod(tsamp*length(kx),gts);
gx = interp1([(0:1:length(gx)-1)*tsamp]',gx,[0:gts:tGrad]','linear','extrap');
gy = interp1([(0:1:length(gy)-1)*tsamp]',gy,[0:gts:tGrad]','linear','extrap');
gz = interp1([(0:1:length(gz)-1)*tsamp]',gz,[0:gts:tGrad]','linear','extrap');

clear kx ky kz
for ll = 1:lshots
    Kx = cumsum([0; gx(:,ll)])*gts*FOVref*gambar;
    Ky = cumsum([0; gy(:,ll)])*gts*FOVref*gambar;
    Kz = cumsum([0; gz(:,ll)])*gts*FOVref*gambar;
    kxt=interp1([0:gts:gts*(length(Kx)-1)],Kx,[0:tsamp:gts*(length(Kx)-1)])';
    kyt=interp1([0:gts:gts*(length(Kx)-1)],Ky,[0:tsamp:gts*(length(Kx)-1)])';
    kzt=interp1([0:gts:gts*(length(Kx)-1)],Kz,[0:tsamp:gts*(length(Kx)-1)])';
    kx(:,ll)=kxt+kx_pre(:,ll);
    ky(:,ll)=kyt+ky_pre(:,ll);
    kz(:,ll)=kzt+kz_pre(:,ll);
end

%ramp down to zero gradient
gx_post = zeros(lRampDownPts,lshots);
gy_post = zeros(lRampDownPts,lshots);
gz_post = zeros(lRampDownPts,lshots);

for xx = 1:lRampDownPts
    gx_post(xx,:) = gx(end,:)-gx(end,:)*xx/lRampDownPts;
    gy_post(xx,:) = gy(end,:)-gy(end,:)*xx/lRampDownPts;
    gz_post(xx,:) = gz(end,:)-gz(end,:)*xx/lRampDownPts;
end

%combine all the gradients together
gx = [gx_pre; gx; gx_post];
gy = [gy_pre; gy; gy_post];
gz = [gz_pre; gz; gz_post];


%must be between +1 and -1
gx = gx/gAmp;
gy = gy/gAmp;
gz = gz/gAmp;

%invert Z gradient to match interpretation in sequence
gz = -gz;

lgradpts = size(gx,1);		
%Set FOV limits to prevent exceeding gradient amplitude and slew limits

%check to make sure gradient limits are not exceeded TODO
FOVminX = floor(FOVref/gAmp);
FOVminY = floor(FOVref/gAmp);
FOVminZ = floor(FOVref/gAmp);

%check to make sure slew rate is not exceeded TODO

gInfo = [gAmp; FOVref; ladc; tsamp; FOVminX; FOVminY; FOVminZ; lshots; Nx; Ny; Nz; tCenter; lgradpts; lRampUpPts; lRampDownPts];

fidx = fopen('gradpulse.gx','w');
fidy = fopen('gradpulse.gy','w');
fidz = fopen('gradpulse.gz','w');
fidinfo = fopen('gradpulse.info','w');

fwrite(fidx,gx(:),'float');
fwrite(fidy,gy(:),'float');
fwrite(fidz,gz(:),'float');
fwrite(fidinfo,gInfo,'float');
fclose('all')

save kspace_traj kx ky kz ladc tsamp lshots Nx Ny Nz tCenter

end
