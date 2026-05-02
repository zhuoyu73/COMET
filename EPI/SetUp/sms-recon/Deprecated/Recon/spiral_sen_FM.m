% spiral_recon_wrapper


TElist = [0,1,2,3]*1e-3;
%TElist = [0,1,10,15]*1e-3;  % Echo Times for Cooling

fname = dir('meas*')
fname = fname(1).name;

sprintf('Creating FM and SENSE map from %s. \n',fname)

curdir = pwd;


recoInfo = setup_recon(fname)
spiral_recon_run(recoInfo);
cd(curdir)

sen = create_sen_map('.', recoInfo);
sprintf('SENSE map created from time 1. \n')

recoInfo.rec_type = 7;
spiral_recon_run(recoInfo);


cd(curdir)
FM = create_field_map_telist('.',recoInfo,TElist);
sprintf('Field map created using TElist: \n')
TElist

