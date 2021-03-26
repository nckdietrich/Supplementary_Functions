function zg_new = interpolate_ZG(zg_old, ilev, lev)
%% Interpolate ZG
%%% Exponentially interpolate ZG values from ilev to lev
%%% All three inputs are column vectors

fcn = @(b, ilev) b(1).*exp(b(2).*ilev);
zg_new = length(zg_old);
for i = 2:length(zg_old)
    zg_part = zg_old(i-1:i);
    ilev_part = ilev(i-1:i);
    lev_part = lev(i-1:i);
    
    b1 = log(zg_part(2)/zg_part(1))/(ilev_part(2) - ilev_part(1));
    a1 = zg_part(1)*exp(-b1*ilev_part(1));
    B = [a1; b1];
    
    zg_est_part_lev = fcn(B, lev_part);
    if i == 2
        zg_new(1:2) = zg_est_part_lev;
    else
        zg_new(i) = zg_est_part_lev(2);
    end
end
end