%% MPC Simulation Script
clear;
close all;
addpath("data\");
load("20-obs-swarm.mat");
load("20-obs.mat")
openfig("20-obs-map.fig");

%% Load Parameters
CarModel = 'Drone';
CarNumber = length(DroneSwarm);
% CarNumber = 1;

MPC_vars = getMPC_vars(CarModel);
ModelParams=getModelParams(MPC_vars.ModelNo);
% choose optimization interface options: 'Yalmip','CVX','hpipm','quadprog'
MPC_vars.interface = 'quadprog';
nx = ModelParams.nx;% number of states
nu = ModelParams.nu;% number of inputs
N = MPC_vars.N;% prediction horizon
Ts = MPC_vars.Ts;% sanple time

%% Simulation lenght and plotting
ts = DroneSwarm(1).ts;
simN = ts(end)*10;
%0=no plots, 1=plot predictions
plotOn = 1;
%0=real time iteration, 1=fixed number of QP iterations, 2=fixed number of damped QP iterations
QP_iter = 0;
%% Fit spline to track
% TODO spline function only works with regular spaced points.
% Fix add function which given any center line and bound generates equlally
% space tracks.
for k=1:CarNumber 
    num = 100;
    sampleTime = linspace(ts(1),ts(end),num);
    track.center = ppval(DroneSwarm(k).pp,sampleTime);
    % plot3(track.center(1,:), track.center(2,:),  track.center(3,:), '--r','LineWidth',1);
    
    track.pp = DroneSwarm(k).pp;
    track.dpp = fnder(DroneSwarm(k).pp,1);
    track.ddpp = fnder(DroneSwarm(k).pp,2);


    % store all data in one struct
    Drone(k).TrackMPC = track;
    Drone(k).ts = ts;

    %% Define starting position
    x0 = zeros(nx,1);
    x0(1:3) = ppval(Drone(k).TrackMPC.pp,ts(1));% current position
    x0(4:6) = ppval(Drone(k).TrackMPC.dpp,ts(1));% current velocity
    Drone(k).x0=x0;
    x = repmat(x0,1,N+1); % all points identical to current measurment
    % first inital guess, all points on centerline aligned with centerline
    for i = 1:N+1
        time_next = (i-1)*Ts;
        x(:,i) = [ppval(Drone(k).TrackMPC.pp,time_next);...
            ppval(Drone(k).TrackMPC.dpp,time_next)]; %driving straight with vx0, and correct theta progress
    end
    u = zeros(nu,N); % zero inputs
    uprev = zeros(nu,1); % last input is zero
    Drone(k).x = x;
    Drone(k).u = u;
    Drone(k).uprev = uprev;
    %% initializtion
    % solve problem 5 times without applying input
    % inspiered by sequential quadratic programming (SQP)
    for i = 1:5
        % formulate MPCC problem and solve it
        Iter_damping = 0; % 0 no damping
        [x_up, u_up, exitflag,info] = optimizer_mpc(0,Drone(k).TrackMPC,...,
            MPC_vars,ModelParams,x,u,x0,uprev,k,Drone);
        x = Iter_damping*x + (1-Iter_damping)*x_up;
        u = Iter_damping*u + (1-Iter_damping)*u_up;
        Drone(k).x = x;
        Drone(k).u = u;
        
    end
end
%% Simulation
mindis_log = [];
maxdis_log = [];
meandis_log = [];
err_log = [];
count = 0;
for k1=1:CarNumber
    Drone(k).flag = 0;
end
for i = 1: simN*10
    count = count + 1;
    for k = 1:CarNumber
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% MPCC-Call %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % augment state and inputs by shifting previus optimal solution
        [Drone(k).x,Drone(k).u] = augState(Drone(k).x,Drone(k).u,Drone(k).x0,MPC_vars,ModelParams);
        %  formulate MPCC problem and solve it
        if QP_iter == 0
            Iter_damping = 0;
            [x_up, u_up, exitflag,info] = optimizer_mpc(i*MPC_vars.Ts,Drone(k).TrackMPC,...,
            MPC_vars,ModelParams,Drone(k).x, Drone(k).u, Drone(k).x0, Drone(k).uprev,k,Drone);
            Drone(k).x = Iter_damping*Drone(k).x + (1-Iter_damping)*x_up;
            Drone(k).u = Iter_damping*Drone(k).u + (1-Iter_damping)*u_up;
        elseif QP_iter == 1
            % doing multiple "SQP" steps
            for k1 = 1:2
                [x_up, u_up, exitflag,info] = optimizer_mpcc(Drone(k).TrackMPC,MPC_vars,ModelParams, Drone(k).x, Drone(k).u, Drone(k).x0, Drone(k).uprev,k,Drone);
                 Drone(k).x = Iter_damping*Drone(k).x + (1-Iter_damping)*x_up;
                Drone(k).u = Iter_damping*Drone(k).u + (1-Iter_damping)*u_up;
            end
        elseif QP_iter == 2
            % doing multiple damped "SQP" steps
            for k1 = 1:2
                Iter_damping = 0.75; % 0 no damping
                [x_up, u_up, exitflag,info] = optimizer_mpcc(Drone(k).TrackMPC,MPC_vars,ModelParams, Drone(k).x, Drone(k).u, Drone(k).x0, Drone(k).uprev,k,Drone);
                Drone(k).x = Iter_damping*Drone(k).x + (1-Iter_damping)*x_up;
                Drone(k).u = Iter_damping*Drone(k).u + (1-Iter_damping)*u_up;
            end
        else
            error('invalid QP_iter value')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% simulate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Drone(k).x0 = SimTimeStep(Drone(k).x(:,1),Drone(k).u(:,1),Ts,ModelParams)';
        Drone(k).uprev = Drone(k).u(:,1);

        if norm(Drone(k).x0(1:3) - ppval(Drone(k).TrackMPC.pp,ts(end))) < 0.5
            Drone(k).flag = 1;
        end
       
    end

    sumFlag = 0;
    for k1=1:CarNumber
        sumFlag = sumFlag + Drone(k1).flag;
    end
    if sumFlag >= CarNumber
       break;
    end

    
    mindis_tmp = 1000;
    maxdis_tmp = 0;
    meandis_tmp = 0;
    count_k1 = 0;
    for k1 = 1:CarNumber
        pos1 = Drone(k1).x0(1:3);
        for k2=k1+1:CarNumber
            pos2 = Drone(k2).x0(1:3);
            dis = norm(pos1 - pos2);
            if dis < mindis_tmp
                mindis_tmp = dis;
            end
            if dis > maxdis_tmp
                maxdis_tmp = dis;
            end
            meandis_tmp = meandis_tmp + dis;
            count_k1 = count_k1 + 1;
        end
        pose_des = ppval(Drone(k1).TrackMPC.pp,i*MPC_vars.Ts);
        pose_err = norm(pose_des - pos1);
        err_log(k1,count) = pose_err;
    end
    meandis_tmp = meandis_tmp / count_k1;
    mindis_log = [mindis_log mindis_tmp];
    maxdis_log = [maxdis_log maxdis_tmp];
    meandis_log = [meandis_log meandis_tmp];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% plotting and logging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotOn == 1
        % figure(2), 
        for k = 1:CarNumber
            Drone(k).pose(:,count) = Drone(k).x0(:,1);
            Drone(k).speed(count) = norm(Drone(k).x0(4:6,1));
            % scatter3(Drone(k).x(1,1),Drone(k).x(2,1),Drone(k).x(3,1),'filled','or');
            % hold on;
            % plot3(Drone(k).x(1,:),Drone(k).x(2,:),Drone(k).x(3,:),'-k','LineWidth',1);
            track.center = ppval(DroneSwarm(k).pp,sampleTime);
            % plot3(track.center(1,:), track.center(2,:),  track.center(3,:), '-r','LineWidth',1);
        end
        axis equal;
        view([-1.678312494634496e+02,60.86259206553521])
        pause(0.01);
        hold off;
    end


   
end

figure(1), hold on; axis equal;
% show(map3D);
% colormap gray;
for k=1:CarNumber
    max_speed = max(Drone(k).speed);
    for k1=1:length(Drone(k).speed)
        RGB(k1,1:3) = map_to_rgb(Drone(k).speed(k1)/max_speed);
    end
    scatter3(Drone(k).pose(1,:), Drone(k).pose(2,:), Drone(k).pose(3,:),2,...
        RGB);
end
title('');
xlabel("x (m)", 'FontSize',8);
ylabel("y (m)", "FontSize",8);
zlabel("z (m)", 'FontSize',8);
xlim([0 240.0])
ylim([0 200.0])
zlim([0.00 40.0])
mean_tmp = 0;
for k=1:CarNumber
    mean_tmp = mean_tmp + mean(Drone(k).speed);
end
mean_speed = mean_tmp / CarNumber;