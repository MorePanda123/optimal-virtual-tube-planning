function polys = minimum_snap_single_axis_simple(waypts,ts,n_order,v0,a0,ve,ae,n_obj,path)
% input:
% waypts:一维位置序列
% ts:时间分配序列
% n_order:多项式阶数
% v0:初始速度
% a0:初始加速度
% ve:最终速度
% ae:最终加速度
% n_obj:目标函数优化阶数
% path:Nx2 N表示段数，前1列表示位置，第2列表示半径
% output:
% polys:m*n矩阵，其中m表示多项式系数数量，n表示多项式分段数

p0 = waypts(1);%初始位置
pe = waypts(end);%终点位置

n_poly = length(waypts)-1;%划分多项式段数
n_coef = n_order+1;%多项式系数数量为点数+1

% compute Q 目标函数
Q_all = [];
for i=1:n_poly
    %blkdiag:分块对角矩阵
    Q_all = blkdiag(Q_all,computeQ(n_order,n_obj,ts(i),ts(i+1)));
end
b_all = zeros(size(Q_all,1),1);

Aeq = zeros(4*n_poly+2,n_coef*n_poly);
beq = zeros(4*n_poly+2,1);

% start/terminal pva constraints  (6 equations)
% 初始端、末端的时间、加速度、速度、位置的约束
Aeq(1:3,1:n_coef) = [calc_tvec(ts(1),n_order,0);
                     calc_tvec(ts(1),n_order,1);
                     calc_tvec(ts(1),n_order,2)];
Aeq(4:6,n_coef*(n_poly-1)+1:n_coef*n_poly) = ...
                    [calc_tvec(ts(end),n_order,0);
                     calc_tvec(ts(end),n_order,1);
                     calc_tvec(ts(end),n_order,2)];
beq(1:6,1) = [p0,v0,a0,pe,ve,ae]';

% mid p constraints    (n_ploy-1 equations)
neq = 6;
for i=1:n_poly-1
    neq=neq+1;
    Aeq(neq,n_coef*i+1:n_coef*(i+1)) = calc_tvec(ts(i+1),n_order,0);
    beq(neq) = waypts(i+1);
end

% continuous constraints  ((n_poly-1)*3 equations)
for i=1:n_poly-1
    tvec_p = calc_tvec(ts(i+1),n_order,0);
    tvec_v = calc_tvec(ts(i+1),n_order,1);
    tvec_a = calc_tvec(ts(i+1),n_order,2);
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_p,-tvec_p];
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_v,-tvec_v];
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_a,-tvec_a];
end

% 不等式约束
Aieq = [];
bieq = [];
Aieq = zeros(2*n_poly,n_coef*n_poly);
bieq = zeros(2*n_poly,1);


options = optimoptions('quadprog','MaxIterations',200);
p = quadprog(Q_all,b_all,Aieq,bieq,Aeq,beq,[],[],[],options);

polys = reshape(p,n_coef,n_poly);

end