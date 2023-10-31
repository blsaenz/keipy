function dn = kei_dn(k,start_yr)

    k_len = length(k.day);
    dn = (0:1/24:k_len/24)+k.day(1);
    dn = dn(1:k_len);
    dn = dn + datenum(start_yr,1,1);

end