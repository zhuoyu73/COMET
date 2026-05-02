function ob = removemask( ob )

    ob.is_masked =0;
    ob.mask = [];
    G = ob.G;
    G =removemask(G);
    ob.G =G;
    
end

