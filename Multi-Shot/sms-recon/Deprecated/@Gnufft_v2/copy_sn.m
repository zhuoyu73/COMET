function ob = copy_sn(ob,ob_toCopyFrom)

st = ob.st;
if isobject(ob_toCopyFrom)
    st_tmp = ob_toCopyFrom.st;
    st.sn = st_tmp.sn;
elseif (ob_toCopyFrom ==0)
    st.sn = 0;
end

ob.st = st;
