% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ X,U,dU,info ] = fminconInterface(stage,MPC_vars,ModelParams,Drone)
nx = ModelParams.nx;
nu = ModelParams.nu;
nz = 2*nx+nu;
nxu = nx+nu;
N = MPC_vars.N;
nc = length(stage(1).lg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = zeros(N*nxu+nx,N*nxu+nx);
f = zeros(N*nxu+nx,1);
for i = 1:N+1
    if i <= N
        H_i = blkdiag(stage(i).Qk,stage(i).Rk);
        H((i-1)*nxu+(1:nxu),(i-1)*nxu+(1:nxu)) = H_i;
%         f((i-1)*nxu+(1:nxu)) = stage(i).fk;
    elseif i == N+1
        H_i = stage(i).Qk;
        H((i-1)*nxu+(1:nx),(i-1)*nxu+(1:nx)) = H_i;
%         f((i-1)*nxu+(1:nx)) = stage(i).fk;
    end
    
end

H = 0.5*(H+H');

Aeq = zeros(nx*(N+1),N*nxu+nx);
beq = zeros(nx*(N+1),1);



Aeq(1:nx,1:nx) = eye(nx);
beq(1:nx) = stage(1).x0;


for i = 1:N

Aeq(i*nx+(1:nx),(i-1)*nxu+(1:nz)) = stage(i).Ak;
beq(i*nx+(1:nx)) = stage(i).Bk;
end

A = [];%ones(1,nz*N+nxu);
b = [];%zeros(ng*N,1);
if nc == 0
    A = [];
    b = [];
else
    A = zeros(nc*N,N*nxu+nx);
    b = zeros(nc*N,1);
    for i = 1:N-1
        A((i-1)*nc + (1:nc),i*nxu+(1:nxu)) = stage(i+1).Ck;
        b((i-1)*nc+(1:nc)) = stage(i+1).lg;
    end
end


LB = zeros(nxu*N+nx,1);
UB = zeros(nxu*N+nx,1);

for i=1:N+1
    if i < N+1
        LB((i-1)*nxu+[1:nxu]) = stage(i).lb;
        UB((i-1)*nxu+[1:nxu]) = stage(i).ub;
    else
        LB((i-1)*nxu+(1:nx)) = stage(i).lb(1:nx);
        UB((i-1)*nxu+[1:nx]) = stage(i).ub(1:nx);
    end
   
end

% options = optimoptions(@quadprog,'MaxIterations',100,'Display','off');
% options = optimoptions(@quadprog,'Display','off');
tic
% [z,~,exitflag] = quadprog(H,f,A,b,Aeq,beq,LB,UB,[]);
options = optimoptions(@fmincon,'MaxFunctionEvaluations',100*length(f));
for k=1:N
    x0((k-1)*nxu+1:(k-1)*nxu+nxu,1) = [Drone.x(:,k);Drone.u(:,k)];
end
x0 = [x0;Drone.x(:,end)];
% x0 = 0.5*(LB+UB);
[z,~,exitflag] = fmincon(@(x)objFunc(x,H,f),x0,A,b,Aeq,beq,LB,UB,[],options);
QPtime = toc;

X = zeros(nx,N+1);
U = zeros(nu,N);
dU = [];
if exitflag == 1
    for i = 1:N+1
        X(1:nx,i) = z((i-1)*nxu+[1:nx]);
        if i < N+1
            U(1:nu,i) = z((i-1)*nxu+nx+[1:nu]);
        end
        
    end
end

info.QPtime = QPtime;
if exitflag == 1
    info.exitflag = 0;
else
    info.exitflag = 1;
end

end

function fval = objFunc(x,H,f)
fval = 0.5*x'*H*x + f'*x;
end
