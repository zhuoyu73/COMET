function ob = removemask( ob, mask )

    ob.is_masked =0;
    ob.mask = [];
    
    dims = size(ob);
    if ob.is_masked	
        dims(end) = sum(ob.mask);
    else
        dims(end) = prod(ob.Nd);
    end
    ob.dims = dims;
end

