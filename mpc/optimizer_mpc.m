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

function [x, u, exitflag,info] = optimizer_mpc(time,TrackMPC,MPC_vars,...
    ModelParams,x,u,x0,u0,num,Drone)

    % reshape X such that it can be used in DP
    X = reshape(x,(MPC_vars.N+1)*ModelParams.nx,1);
    % add other cars to X for DP
%     if ~isempty(Y)
%         Y = repmat(Y,MPC_vars.N+1,1);
%     end
%     X = [X,Y];
    % get boarder points using DP based algorithm
%     [new_b,~] = getNewBorders(TrackMPC.traj,TrackMPC.borders,TrackMPC.track_center,X,n_cars,1,MPC_vars,ModelParams);
    % formulate MPCC problem and solve it given the DP bundaries
    [xpred, upred,dupred,info] = mao_getMPCmatrices(time,MPC_vars,...
        ModelParams,x,u,x0,u0, num, Drone);

    if (info.exitflag == 0)
        exitflag = 0;
        x = xpred;
        u = upred;
    else
        exitflag = 1;
%         x(:,1:end-1) = x(:,2:end);
%         u(:,1:end-1) = u(:,2:end);
        % x and u stay identical to the initial guess
    end
end