% spiral_recon_wrapper

dirnm_FM1 = '4ms'
dirnm_FM2 = '4p5ms'

curdir = pwd;

if 1
  cd(dirnm_FM1)
  recoInfoFM1 = setup_recon;
  spiral_recon_run(recoInfoFM1);
  cd(curdir)

  sen = create_sen_map(dirnm_FM1,recoInfoFM1);

  cd(curdir)
  cd(dirnm_FM1)
  recoInfoFM1 = setup_recon;
  recoInfoFM1.rec_type = 7;
  spiral_recon_run(recoInfoFM1);

  cd(curdir)
  cd(dirnm_FM2)
  recoInfoFM2 = setup_recon;
  recoInfoFM2.rec_type = 7;
  spiral_recon_run(recoInfoFM2);

  cd(curdir)
  FM = create_field_map(dirnm_FM1,dirnm_FM2,recoInfoFM1,recoInfoFM2);

end
if 0
% recoInfo = setup_recon;
recoInfo = setup_recon([],[],[],[],[],3);  
recoInfo.rec_type = 8;
  recoInfo.shots_to_use = recoInfo.nl;
  %spiral_recon_run(recoInfo);
 spiral_saveraw(recoInfo);


%  !mv Recon ReconSENSE1

%  recoInfo.shots_to_use = 1;
%  spiral_recon_run(recoInfo);
%  !mv Recon ReconSENSE3

load recoInfo
cd Recon
fnames = dir('*.mat');

for ii = 1:length(fnames)
    eval(sprintf('load %s',fnames(ii).name))
    name_to_use = fnames(ii).name;
    name_to_use = name_to_use(1:end-4);
    write_ana(name_to_use,permute(img,[2 1 3]),recoInfo.szx, recoInfo.szy, recoInfo.szz,'FLOAT',1);
end





end
