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

function [Ad, Bd, gd]=DiscretizedLinearizedModel(Xbar_k,Ubar_k,ModelParams,Ts)
% returns the discretized, linearized model about (Xbar_k,Ubar_k)
% s.t. x(k+1) = Ad*x(k) + Bd*u(k) + gd
%
% Ad: [nx]x[nx], Bd: [nx]x[nu], gd: [nx]x[1]

  
if ModelParams.ModelNo==1 || ModelParams.ModelNo==3
    l = ModelParams.l;
    Ad = [1  0  0 Ts         0            0     0;
        0  1  0  0         Ts            0     0;
        0  0  1  0          0           Ts     0;
        0  0  0 1-Ts*l  0            0      0;
        0  0  0 0       1-Ts*l       0      0;
        0  0  0 0          0       1-Ts*l   0;
        0  0  0 0          0            0       1];
    Bd = [0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        Ts*l 0 0 0;
        0 Ts*l 0 0;
        0 0 Ts*l 0;
        0 0 0 Ts];
    gd = zeros(ModelParams.nx,1);
elseif ModelParams.ModelNo == 4
    Ad = [1 0 0 Ts 0 0;
          0 1 0 0 Ts 0;
          0 0 1 0 0 Ts;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1];
    Bd = [0 0 0;
        0 0 0;
        0 0 0;
        Ts 0 0;
        0 Ts 0;
        0 0 Ts];
    gd = [];
elseif ModelParams.ModelNo == 2
    Ad = [
        1 0 0 Ts 0 0 0 0 0 -1 0 0 0 0 0;
        0 1 0 0 Ts 0 0 0 0 0 -1 0 0 0 0;
        0 0 1 0 0 Ts 0 0 0 0 0 -1 0 0 0;
        0 0 0 1 0 0 Ts 0 0 0 0 0 -1 0 0;
        0 0 0 0 1 0 0 Ts 0 0 0 0 0 -1 0;
        0 0 0 0 0 1 0 0 Ts 0 0 0 0 0 -1];
    Bd = [];
    gd = [];
end
    



% in fact, the above can be done using only the physical states/inputs, 
% then can add    Ad(end+1,end+1)=1; Bd(end+1,end+1)=Ts; 
    
end