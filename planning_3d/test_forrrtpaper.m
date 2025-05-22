clear;close all;
addpath("tubeRRTStar\");
addpath("data\");
%% open map figure
openfig("60-obs-map.fig");
load("60-obs.mat");
%% generate paths and optimize trajectories
startPoint = [10 110 15 15];
sp1 = [10 100 15];
sp2 = [10 120 15];
sp3 = [10 110 5];

goalPoint = [220 100 15 15];
gp1 = [220 90 15];
gp2 = [220 110 15];
gp3 = [220 100 5];
% path points
path = path;
path(1,:) = [];
% path(end-1,:) = [220 100 15 15];
path(end,:) = [];
tubePaths = genTubePath2(path);
wp1 = [sp1;real(tubePaths(1).path);gp1];
wp2 = [sp2;real(tubePaths(2).path);gp2];
wp3 = [sp3;real(tubePaths(3).path);gp3];

hold on;
% plot3(path(:,1),path(:,2),path(:,3));
plot3(wp1(:,1),wp1(:,2),wp1(:,3));
plot3(wp2(:,1),wp2(:,2),wp2(:,3));
plot3(wp3(:,1),wp3(:,2),wp3(:,3));

% optimal trajectories
numSamples = 100;
v_mean = 4;

timePoint1 = 0;
timePoint2 = 0;
timePoint3 = 0;
timePoint4 = 0;
for k=1:length(wp1(:,1))-1
    t_tmp = norm(wp1(k+1,:) - wp1(k,:))/v_mean;
    timePoint1 = [timePoint1 timePoint1(end)+t_tmp];
    t_tmp = norm(wp2(k+1,:) - wp2(k,:))/v_mean;
    timePoint2 = [timePoint2 timePoint2(end)+t_tmp];
    t_tmp = norm(wp3(k+1,:) - wp3(k,:))/v_mean;
    timePoint3 = [timePoint3 timePoint3(end)+t_tmp];
 
end
timePoints = (timePoint3 + timePoint2 + timePoint1) / 3;
timePoints = timePoints / (timePoints(end) - timePoints(1)) * 90;

[traj1,~,~,~,pp1,~,~] = minjerkpolytraj(wp1',timePoints,numSamples);
[traj2,~,~,~,pp2,~,~] = minjerkpolytraj(wp2',timePoints,numSamples);
[traj3,~,~,~,pp3,~,~] = minjerkpolytraj(wp3',timePoints,numSamples);

plot3(traj1(1,:),traj1(2,:),traj1(3,:),'LineWidth',2);
plot3(traj2(1,:),traj2(2,:),traj2(3,:),'LineWidth',2);
plot3(traj3(1,:),traj3(2,:),traj3(3,:),'LineWidth',2);

num = 100;
sampleTime = linspace(timePoints(1),timePoints(end),num);
count = 1;
for k1=0.1:0.2:1-0.1
    for k2=0.1:0.2:1-k1
        k3 = 1 - k1 - k2;
        DroneSwarm(count).ts = timePoints;
        DroneSwarm(count).pp = pp1;
        DroneSwarm(count).pp.coefs = k1*pp1.coefs + k2*pp2.coefs + k3*pp3.coefs;
        
        curve_tmp = ppval(DroneSwarm(count).pp,sampleTime);
        plot3(curve_tmp(1,:),curve_tmp(2,:),curve_tmp(3,:),'LineWidth',2,'Color','g');
        count = count + 1;
    end
end
