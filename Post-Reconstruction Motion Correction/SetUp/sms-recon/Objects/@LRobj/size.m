function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class
length_total = 0;
list=ob.A_index_list;
for kk = 1:ob.ntp
	aa=list(kk);
    sizeA = size(ob.A{aa});
    length_total = length_total + sizeA(1);
end
dim = [length_total,sizeA(2)*ob.R] ;


