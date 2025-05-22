function [flag,path_inserted] = PathInsert(path,gen_curve,num)
global bw1 r_max
flag = 0;
for k = 1:length(gen_curve.origin)
    count = ceil(k/num);
    x_b = round(gen_curve.origin(1,k));
    y_b = round(gen_curve.origin(2,k));
    if bw1(y_b,x_b) == 1
        path_ob = [x_b;y_b];
        vec_0 = path_ob - path(1:2,count);
        vec_1 = path(1:2,count+1) - path(1:2,count);
        vec_2 = (vec_0' * vec_1) / norm(vec_1)^2 * vec_1 + path(1:2,count);
        path_tmp = [vec_2;r_max];
        path_inserted = [path(:,1:count) path_tmp path(:,count+1:end)];
        flag = 1;
        break;
    else
        flag =0;
    end
end
if flag == 0 
    path_inserted = path;
end
end
