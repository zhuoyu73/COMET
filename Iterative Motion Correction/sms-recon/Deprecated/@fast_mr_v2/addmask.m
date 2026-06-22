function ob = addmask( ob, mask )

    ob.is_masked =1;
    ob.mask = mask;
    G = ob.G;
    G =addmask(G,mask);
    ob.G =G;
end

