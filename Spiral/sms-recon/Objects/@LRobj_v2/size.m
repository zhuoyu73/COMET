function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class
length_total = 0;
for kk = 1:ob.ntp
    sizeA = size(ob.A{kk});
    length_total = length_total + sizeA(1);
end
dim = [length_total,sizeA(2)*ob.R] ;


