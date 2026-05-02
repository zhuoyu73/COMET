function ob = copy_sn(ob,ob_toCopyFrom)

if isobject(ob_toCopyFrom)
    ob.G = copy_sn(ob.G,ob_toCopyFrom.G);

elseif (ob_toCopyFrom ==0)
    ob.G = copy_sn(ob.G,0);
end


