function dim = size(ob)

nobj = ob.n;
A = ob.A{1};
s_map = ob.smap;
sA = size(A);
dim = [nobj*sA(1)*size(s_map, 2) sA(2)];
