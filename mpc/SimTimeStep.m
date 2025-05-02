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

function xp=SimTimeStep(x,u,Ts,ModelParams)
%x state
%u input
%Ts sampling time
x0=x;

[~,inivt]=ode45(@(t,x)fx_bicycle(t,x,u,ModelParams),[0 Ts],x0);
xp=inivt(end,:);
return

function xdot=fx_bicycle(t,x,u,ModelParams)

if ModelParams.ModelNo == 2
    l = 2;
    xdot = [x(4);
        x(5);
        x(6);
        -l*(x(4)-u(1));
        -l*(x(5)-u(2));
        -l*(x(6)-u(3))];
elseif ModelParams.ModelNo == 1
    xdot = [x(4);
        x(5);
        x(6);
        u(1);
        u(2);
        u(3)];
end
    
return