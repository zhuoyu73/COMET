function [kspace, grads_mTm] = gen_nav_SpiralInOut_readout(N,multibandFactor, Rxy,Rz, FOV,FOVz, gAmp, gSlew, tSamp)

    
% We design the spiral in case by default.    
[kspace, grads_mTm] = gen_nav_readout(N, multibandFactor, Rxy*2, Rz, FOV, FOVz, gAmp, gSlew, tSamp);

grads_mTm_Rev = grads_mTm(end:-1:1,:);
grads_mTm_Rev(:,3) = -1*grads_mTm(end:-1:1,3); % Deal with Z Separately

grads_mTm = vertcat(grads_mTm,grads_mTm_Rev);

kspace_rev(:,1) = -1*kspace(end:-1:1,1);
kspace_rev(:,2) = -1*kspace(end:-1:1,2);
kspace_rev(:,3) =    kspace(end:-1:1,3);

kspace = vertcat(kspace, kspace_rev);
centerOfKSpace = length(grads_mTm_Rev) + 1;
filename = ['Rxy' num2str(Rxy) 'Rz' num2str(Rz) 'Nxy' num2str(N) 'Nz' num2str(multibandFactor) '_InOut.xml'];
writeExternalGradFile(filename, grads_mTm, FOV, FOVz, centerOfKSpace,N);
save(['Rxy' num2str(Rxy) 'Rz' num2str(Rz) 'Nxy' num2str(N) 'Nz' num2str(multibandFactor) '_InOut.mat'],'kspace','grads_mTm');
end
