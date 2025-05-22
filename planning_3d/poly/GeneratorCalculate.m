function [gen_curve] = GeneratorCalculate(polys_x, polys_y, polys_z, ts, num)

len = num * (length(ts)-1);
fx=zeros(1,len);
fx1 = fx;
fx2 = fx;
fy = fx;
fy1 = fx;
fy2 = fx;
fz = fx;
ft = fx;
for i=1:size(polys_x,2)
    tt = linspace(ts(i),ts(i+1),num);
    xx = polys_vals(polys_x,ts,tt,0);
    xx1 =polys_vals(polys_x,ts,tt,1);
    xx2 =polys_vals(polys_x,ts,tt,2);
    yy = polys_vals(polys_y,ts,tt,0);
    yy1 = polys_vals(polys_y,ts,tt,1);
    yy2 = polys_vals(polys_y,ts,tt,2);
    zz = polys_vals(polys_z,ts,tt,0);
    fx(num*(i-1)+1:num*i) = xx;
    fx1(num*(i-1)+1:num*i) = xx1;
    fx2(num*(i-1)+1:num*i) = xx2;
    fy(num*(i-1)+1:num*i) = yy;
    fy1(num*(i-1)+1:num*i) = yy1;
    fy2(num*(i-1)+1:num*i) = yy2;
    fz(num*(i-1)+1:num*i) = zz;
    ft(num*(i-1)+1:num*i) = tt;
end

gen_curve.origin = [fx; fy; fz];
gen_curve.d1 = [fx1; fy1];
gen_curve.d2 = [fx2; fy2];
gen_curve.t = ft;

% unit principal normal
% tangent
norm_f1 = sum(gen_curve.d1.^2,1).^0.5;
% tangent = gen_curve.d1./norm_f1;
prin_normal = gen_curve.d2./norm_f1 - sum(gen_curve.d1.* gen_curve.d2,1)./(norm_f1.^3).*gen_curve.d1;
prin_normal = prin_normal./ sum(prin_normal.^2,1).^0.5;
gen_curve.principal_normal = prin_normal;
end