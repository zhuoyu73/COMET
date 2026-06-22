function ob = copy_timeSegmentationParams(ob,ob_toCopyFrom)

if isobject(ob_toCopyFrom)
    ob.we = ob_toCopyFrom.we;
    ob.int = ob_toCopyFrom.int;
elseif (ob_toCopyFrom ==0)
    ob.we = 0;
    ob.int =0;
end


