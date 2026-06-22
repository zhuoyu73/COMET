function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class
dim = [size(ob.C,1)*size(ob.R),ob.N*ob.R] ;


