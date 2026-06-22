function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class
tmp = ob.G;
tmp2 = ob.HH;

ndims1 = size(tmp);
ndims2 = size(tmp2);

dim = [ndims1(1),ndims2(2)];
