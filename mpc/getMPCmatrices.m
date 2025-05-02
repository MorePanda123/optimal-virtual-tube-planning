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

function [X,U,dU,info] = getMPCmatrices(traj,MPC_vars,ModelParams,Xhor,Uhor,x0,u0, num, Drone)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each stage in the horizon compute the necessary system and cost %%%%%
% matricies and safe them in a array of structs 'stage' %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost scaling for numerics
costScale = MPC_vars.costScale;
% init stage struct
stage = struct([]);
%% Generate MPCC problem
% initial state (including previus input for input rate cost)

% state is augmented by the inputs and rate inputs are introduced to
% formulate rate input costs and constraints while retaining a block
% spares formulation
% given x_k+1 = A x_k + B u_k
% do the following state augmentation
% s_k = [x_k,u_k-1], v_k = du_k
% with the following linear system
% s_k+1 = [A B;0 I] s_k + [B;I] v_k

stage(1).x0 = x0;
stage(1).u0 = u0;
if ModelParams.ModelNo == 1
    for i = 1:MPC_vars.N
        Xk = Xhor(:,i);
        Uk = Uhor(:,i);   
        % generate quadratic state(-input) cost
        stage(i).Qk = costScale*generateH(traj,MPC_vars,ModelParams,Xk,i);
        % quadratic rate input cost
        stage(i).Rk = costScale*2*diag([MPC_vars.rdVx; MPC_vars.rdVy; MPC_vars.rdVz; MPC_vars.rdVtheta]);
        % linear state(-input) cost
        stage(i).fk = costScale*generatef(traj,MPC_vars,ModelParams,Xk,i);
        % linearized dynamics
        [stage(i).Ak,stage(i).Bk,stage(i).gk] = getEqualityConstraints(Xk,Uk,MPC_vars,ModelParams);
        % linearized track constraints
        [stage(i).Ck,stage(i).lg] = getInequalityConstraints(MPC_vars,ModelParams,Xk,num,Drone,i);
        % bounds
        [stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);
    end
    % terminal stage
    i = MPC_vars.N+1;
    Xk = Xhor(:,i);
    % generate quadratic state(-input) cost
    stage(i).Qk = costScale*generateH(traj,MPC_vars,ModelParams,Xk,i);
    % quadratic rate input cost
    stage(i).Rk = costScale*2*diag([MPC_vars.rdVx; MPC_vars.rdVy; MPC_vars.rdVz; MPC_vars.rdVtheta]);
    % linear state(-input) cost
    stage(i).fk = costScale*generatef(traj,MPC_vars,ModelParams,Xk,i);
    % linearized track constraints
    [stage(i).Ck, stage(i).lg] = getInequalityConstraints(MPC_vars,ModelParams,Xk,num,Drone,i);
    % bounds
    [stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);
elseif ModelParams.ModelNo == 2
    for i = 1:MPC_vars.N
        Xk = Xhor(:,i);
        Uk = Uhor(:,i);
        % generate quadratic state cost
        stage(i).Qk = diag([MPC_vars.qp,MPC_vars.qp,MPC_vars.qp,...
            MPC_vars.qv,MPC_vars.qv,MPC_vars.qv]);
        % quadratic input cost
        stage(i).Rk = costScale*diag([MPC_vars.qa; MPC_vars.qa;...
            MPC_vars.qa]);
        % linear state(-input) cost
        stage(i).fk = -2*[Xk;Uk]'*blkdiag(stage(i).Qk,stage(i).Rk);
        % linearized dynamics
        [stage(i).Ak,stage(i).Bk,stage(i).gk] = getEqualityConstraints(Xk,Uk,MPC_vars,ModelParams);
        % linearized track constraints
        [stage(i).Ck,stage(i).lg] = getInequalityConstraints(MPC_vars,ModelParams,Xk,num,Drone,i);
        % bounds
        [stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);
    end
    % terminal stage
    i = MPC_vars.N+1;
    Xk = Xhor(:,i);
    % generate quadratic state(-input) cost
    stage(i).Qk = costScale*generateH(traj,MPC_vars,ModelParams,Xk,i);
    % quadratic rate input cost
    stage(i).Rk = costScale*2*diag([MPC_vars.rdVx; MPC_vars.rdVy; MPC_vars.rdVz; MPC_vars.rdVtheta]);
    % linear state(-input) cost
    stage(i).fk = costScale*generatef(traj,MPC_vars,ModelParams,Xk,i);
    % linearized track constraints
    [stage(i).Ck, stage(i).lg] = getInequalityConstraints(MPC_vars,ModelParams,Xk,num,Drone,i);
    % bounds
    [stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);
end
%% Call solver interface
if strcmp(MPC_vars.interface, 'Yalmip')
    % yalmip based interface (very slow)
    [X,U,dU,info] = YalmipInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface, 'CVX')
    % CVX based interface (slow)
    [X,U,dU,info] = CVXInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface, 'hpipm')
    % hpipm interface (prefered)
    [X,U,dU,info] = hpipmInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface, 'quadprog')
    % quadprog interface (replace quadprog with a better solver if possible)
    [X,U,dU,info] = QuadProgInterface(stage,MPC_vars,ModelParams);
else
    error('invalid optimization interface')
end

end

% GENERATING Q
function Qk = generateH(pathinfo,MPC_vars,ModelParams,Xk,i)
    % get linearized contouring and lag errors
    Qtilde = generateQtilde(pathinfo,MPC_vars,ModelParams,Xk,i);
    % add omega regularization
%     if i == MPC_vars.N+1
%         Qtilde(ModelParams.stateindex_omega,ModelParams.stateindex_omega) = MPC_vars.qOmegaNmult*MPC_vars.qOmega;
%     else
%         Qtilde(ModelParams.stateindex_omega,ModelParams.stateindex_omega) = MPC_vars.qOmega;
%     end
    % make Qtilde symetric (not symetric due to numerical issues)
    Qtilde = 0.5 *(Qtilde+Qtilde');
    % Qk = contouring-lag error and real-input cost
    Qk = 2*blkdiag(Qtilde,diag([MPC_vars.rVx,MPC_vars.rVy,MPC_vars.rVz,MPC_vars.rVtheta]));
    % scale cost
    Qk = blkdiag(MPC_vars.invTx,MPC_vars.invTu)*Qk*blkdiag(MPC_vars.invTx,MPC_vars.invTu) + 1e-12*eye(ModelParams.nx+ModelParams.nu);
end

% compute linear contouring and lag errors
function Qtilde = generateQtilde(pathinfo,MPC_vars,ModelParams,Xk,i)
    if i == MPC_vars.N+1
        Q = diag([MPC_vars.qCNmult*MPC_vars.qC, MPC_vars.qC, MPC_vars.qL]);
    else
        Q = diag([MPC_vars.qC, MPC_vars.qC, MPC_vars.qL]);
    end
        
    theta_virt=mod(Xk(end),pathinfo.ppx.breaks(end));
    [grad_eCn, grad_eCb, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams,Xk(1), Xk(2), Xk(3));
    errorgrad = [grad_eCn; grad_eCb; grad_eL]; 
    Qtilde = errorgrad'*Q*errorgrad; 
end

function [grad_eCn, grad_eCb, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams, x_phys, y_phys, z_phys)

    [deCn_dtheta, deCb_dtheta, deL_dtheta, cos_phi_virt, sin_phi_virt,cos_zeta_virt,sin_zeta_virt] = getderror_dtheta(pathinfo, theta_virt, x_phys, y_phys, z_phys);
    
%     grad_eC = [ sin_phi_virt, -cos_phi_virt, zeros(1, ModelParams.nx-3), deC_dtheta];
    grad_eCn = [-sin_zeta_virt*cos_phi_virt,-sin_zeta_virt*sin_phi_virt, cos_zeta_virt, zeros(1, ModelParams.nx-4), deCn_dtheta];
    grad_eCb = [sin_phi_virt, -cos_phi_virt, 0, zeros(1, ModelParams.nx-4), deCb_dtheta];
    grad_eL = [-cos_zeta_virt*cos_phi_virt, -cos_zeta_virt*sin_phi_virt,-sin_zeta_virt, ...
        zeros(1, ModelParams.nx-4), deL_dtheta];
end

function [deCn_dtheta, deCb_dtheta, deL_dtheta, cos_phi_virt, sin_phi_virt,cos_zeta_virt,sin_zeta_virt] = getderror_dtheta(pathinfo, theta_virt, x_phys, y_phys, z_phys)
    dxvirt_dtheta=ppval(pathinfo.dppx,theta_virt); %d x_virt / d theta
    dyvirt_dtheta=ppval(pathinfo.dppy,theta_virt); %d y_virt / d theta
    dzvirt_dtheta=ppval(pathinfo.dppz,theta_virt);
    
    phi_virt=atan2(dyvirt_dtheta,dxvirt_dtheta); %orientation of virtual position
    zeta_virt = atan2(dzvirt_dtheta, sqrt(dxvirt_dtheta^2+dyvirt_dtheta^2));
    
    % virtual positions
    x_virt=ppval(pathinfo.ppx,theta_virt);
    y_virt=ppval(pathinfo.ppy,theta_virt);
    z_virt=ppval(pathinfo.ppz,theta_virt);
    
    % difference in position between virtual and physical
    Dx=x_virt-x_phys;
    Dy=y_virt-y_phys;
    Dz=z_virt-z_phys;

    dphivirt_dtheta=getdphivirt_dtheta(theta_virt,pathinfo);
    dzetavirt_dtheta=getdzetavirt_dtheta(theta_virt,pathinfo);

    cos_phi_virt=cos(phi_virt);
    sin_phi_virt=sin(phi_virt);
    cos_zeta_virt = cos(zeta_virt);
    sin_zeta_virt = sin(zeta_virt);
    
    deL_dtheta = -(sin_zeta_virt * cos_phi_virt * dzetavirt_dtheta +  cos_zeta_virt * sin_phi_virt * dphivirt_dtheta)*Dx...
        + cos_zeta_virt * cos_phi_virt * dxvirt_dtheta...
        - (sin_zeta_virt * sin_phi_virt * dzetavirt_dtheta - cos_zeta_virt * cos_phi_virt * dphivirt_dtheta) * Dy...
        + cos_zeta_virt * sin_phi_virt * dyvirt_dtheta...
        + cos_zeta_virt * dzetavirt_dtheta * Dz...
        + sin_zeta_virt * dzvirt_dtheta;
    
   deCn_dtheta = (dzetavirt_dtheta*cos_zeta_virt*cos_phi_virt - dphivirt_dtheta*sin_zeta_virt*sin_phi_virt)*Dx...
        + dxvirt_dtheta*sin_zeta_virt*cos_phi_virt...
        + (dzetavirt_dtheta*cos_zeta_virt*sin_phi_virt + dphivirt_dtheta*sin_zeta_virt*cos_phi_virt)*Dy...
        + dyvirt_dtheta*sin_zeta_virt*sin_phi_virt...
        + dzetavirt_dtheta*sin_zeta_virt*Dz - dzvirt_dtheta*cos_zeta_virt;
    
    deCb_dtheta = - dphivirt_dtheta*cos_phi_virt*Dx - dxvirt_dtheta*sin_phi_virt...
        - dphivirt_dtheta*sin_phi_virt*Dy...
        +dyvirt_dtheta*cos_phi_virt;
     
end

function dphivirt_dtheta=getdphivirt_dtheta(theta_virt,pathinfo)
    % computes {d phi_virt / d theta} evaluated at theta_k

    dxdth=ppval(pathinfo.dppx,theta_virt); %d x_virt / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); %d y_virt / d theta
    d2xdth2=ppval(pathinfo.ddppx,theta_virt); %d2 x_virt / d theta2
    d2ydth2=ppval(pathinfo.ddppy,theta_virt); %d2 y_virt / d theta2

    numer=dxdth*d2ydth2 - dydth*d2xdth2;
    denom=dxdth^2 + dydth^2;

    dphivirt_dtheta=numer/denom;
end

function dzetavirt_dtheta=getdzetavirt_dtheta(theta_virt,pathinfo)
    dxdth=ppval(pathinfo.dppx,theta_virt); %d x_virt / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); %d y_virt / d theta
    dzdth=ppval(pathinfo.dppz,theta_virt); 
    
    d2xdth2=ppval(pathinfo.ddppx,theta_virt); %d2 x_virt / d theta2
    d2ydth2=ppval(pathinfo.ddppy,theta_virt); %d2 y_virt / d theta2
    d2zdth2=ppval(pathinfo.ddppz,theta_virt); 
    
    numer=d2zdth2*(dxdth^2+dydth^2) - dzdth*(dxdth*d2xdth2+dydth*d2ydth2);
    denom=(dxdth^2+dydth^2+dzdth^2)*sqrt(dxdth^2+dydth^2);
    
    dzetavirt_dtheta = numer / denom;
end

% GENERATING f
function f = generatef(pathinfo,MPC_vars,ModelParams,Xk,i)

    x_phys = Xk(1);
    y_phys = Xk(2);
    z_phys = Xk(3);

    theta_virt=mod(Xk(end),pathinfo.ppx.breaks(end));
    [eCn, eCb, eL] = getErrors(pathinfo, theta_virt,x_phys,y_phys, z_phys);
    e=[eCn; eCb; eL];
    [grad_eCn, grad_eCb, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams, x_phys, y_phys, z_phys);
    grad_e = [grad_eCn; grad_eCb; grad_eL];
    
    if i == MPC_vars.N+1
        Q = diag([MPC_vars.qC,MPC_vars.qC, MPC_vars.qL]);
    else
        Q = diag([MPC_vars.qC,MPC_vars.qC, MPC_vars.qL]);
    end
  
    fx=2*e'*Q*grad_e - 2*Xk'*grad_e'*Q*grad_e;
    fT = [fx, zeros(1,ModelParams.nu-1), -MPC_vars.qVtheta];
    f=fT';
    
    f = blkdiag(MPC_vars.invTx,MPC_vars.invTu)*f;
end

function [eCn, eCb, eL] = getErrors(pathinfo, theta_virt,x_phys,y_phys, z_phys)
    dxdth=ppval(pathinfo.dppx,theta_virt); % d x / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); % d y / d theta
    dzdth=ppval(pathinfo.dppz,theta_virt); 

    % virtual positions
    x_virt=ppval(pathinfo.ppx,theta_virt);
    y_virt=ppval(pathinfo.ppy,theta_virt);
    z_virt=ppval(pathinfo.ppz,theta_virt);
    
    Dx = x_virt - x_phys;
    Dy = y_virt - y_phys;
    Dz = z_virt - z_phys;
    
    phi_virt=atan2(dydth,dxdth);
    zeta_virt = atan2(dzdth, sqrt(dxdth^2+dydth^2));
    
    % define these to reduce calls to trig functions
    sin_phi_virt = sin(phi_virt);
    cos_phi_virt = cos(phi_virt);
    cos_zeta_virt = cos(zeta_virt);
    sin_zeta_virt = sin(zeta_virt);

    % contouring and lag error estimates
%     eC = -sin_phi_virt*(x_virt - x_phys) + cos_phi_virt*(y_virt - y_phys);
   eCn = sin_zeta_virt*cos_phi_virt*Dx + sin_zeta_virt*sin_phi_virt*Dy - cos_zeta_virt*Dz;
   eCb = -sin_phi_virt*Dx + cos_phi_virt*Dy;
%     eL =  cos_phi_virt*(x_virt - x_phys) + sin_phi_virt*(y_virt - y_phys);
    eL = cos_zeta_virt * cos_phi_virt * (x_virt - x_phys) + cos_zeta_virt * sin_phi_virt * (y_virt - y_phys)...
        + sin_zeta_virt * (z_virt - z_phys);
   
end

% EQUALITY CONSTRAINTS
function [Ak,Bk,gk] = getEqualityConstraints(Xk,Uk,MPC_vars,ModelParams) 

    nx = ModelParams.nx;
    nu = ModelParams.nu;
    % linearize and discretize nonlinear bicycle model
    [Ad, Bd, gd]=DiscretizedLinearizedModel(Xk,Uk,ModelParams,MPC_vars.Ts);
    % constructing augmented system with state-input scaling
    Ak = [MPC_vars.Tx*Ad*MPC_vars.invTx MPC_vars.Tx*Bd*MPC_vars.invTu; zeros(nu,nx) eye(nu)];
    Bk = [MPC_vars.Tx*Bd*MPC_vars.invTu;eye(nu)];
    gk = [MPC_vars.Tx*gd;zeros(nu,1)];
    
end

% INEQUALITY CONSTRAINTS
function [Ck, lg] = getInequalityConstraints(MPC_vars,ModelParams,Xk,num,Drone,N_step)
    nx = ModelParams.nx;
    nu = ModelParams.nu;
    rho = ModelParams.radius;
    if num == 1
        Ck=[];
        lg = [];
    else
        Ck=zeros(num-1,nx+nu);
        lg = zeros(num-1,1);
        
        for k =1:num-1
            p_o = Drone(k).x(1:3,N_step);
            
            c_xy = (Xk(1:3) - p_o)/norm(Xk(1:3) - p_o) * rho;
            Ck(k,1:3) = c_xy';
            lg(k) =  2*rho^2 + c_xy'*p_o;
        end
    end
    
    
end

% BOUNDS
function [lb, ub]=getBounds(MPC_vars,ModelParams)

lb = MPC_vars.bounds(:,1);
ub = MPC_vars.bounds(:,2);


end
