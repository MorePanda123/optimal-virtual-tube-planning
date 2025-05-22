clear;close all;
addpath("tubeRRTStar\");
addpath("poly\");
%% generate 3-D occupancy map
map3D = occupancyMap3D(0.1);

pose = [ 0 0 0 1 0 0 0];
randNum = 40;
maxRange = 200;
pointsx = rand(randNum,1)*maxRange+20;
pointsy = rand(randNum,1)*maxRange;
points1 = [pointsx pointsy, ones(randNum,1)];
points2 = [pointsx pointsy, 11*ones(randNum,1)];
points3 = [pointsx pointsy, 21*ones(randNum,1)];
points = [points1;points2;points3];
pointsFloor = [];
for x = 0:1:300
    for y = 0:1:300
        pointsFloor = [pointsFloor; x, y, -0.5];
    end
end

insertPointCloud(map3D,pose,points,maxRange*1.5);
insertPointCloud(map3D,pose,pointsFloor,maxRange*1.5);
map3D.FreeThreshold = map3D.OccupiedThreshold;

xLimMin = 0;
xLimMax = 240;
yLimMin = 0;
yLimMax = 200;
zLimMin = 0;
zLimMax = 30;
%% path planning
startPose = [10 110 7];  % [x y z ]
goalPose = [220 100 15];

% ----------------RRT*-----------------
setting.useIntVol = false;
setting.ContinueAfterGoalReached = true;
setting.MaxConnectionDistance = 10;
setting.GoalBias = 15;
setting.MaxNumTreeNodes = 1000;
setting.xLim = [xLimMin xLimMax];
setting.yLim = [yLimMin yLimMax];
setting.zLim = [zLimMin zLimMax];
setting.useIntVol = true;
setting.costIntVol = 0.015;
[path,tree] = planTubeRRTStar(startPose, goalPose, map3D, setting);
%% plot
figure,
hold on; axis equal;
show(map3D)
colormap gray
xlim([0 240.0])
ylim([0 200.0])
zlim([0.00 60.0])
scatter3(startPose(1),startPose(2),startPose(3),"g",'filled');
scatter3(goalPose(1),goalPose(2),goalPose(3),"g",'filled');

plot3(path(:,1),path(:,2),path(:,3),'-g','linewidth',2);
for k=1:length(path(:,1))
   showCorridor(path(k,4), path(k,1:3)); 
end

 tubePaths = genTubePath2(path);
radius = path(:,4);
%% optimize curve
waypoints1 = tubePaths(1).path;
nwaypoints = length(waypoints1(:,1));
distance = 0;
for k=2:nwaypoints
    distance = distance + norm(waypoints1(k,:) - waypoints1(k-1,:));
end
uavSpeed = 10;
total_time = distance / uavSpeed;
ts1 = arrangeT(waypoints1',total_time);

waypoints2 = tubePaths(2).path;
nwaypoints = length(waypoints1(:,2));
distance = 0;
for k=2:nwaypoints
    distance = distance + norm(waypoints2(k,:) - waypoints2(k-1,:));
end
uavSpeed = 10;
total_time = distance / uavSpeed;
ts2 = arrangeT(waypoints2',total_time);

waypoints3 = tubePaths(3).path;
nwaypoints = length(waypoints3(:,2));
distance = 0;
for k=2:nwaypoints
    distance = distance + norm(waypoints3(k,:) - waypoints3(k-1,:));
end
uavSpeed = 10;
total_time = distance / uavSpeed;
ts3 = arrangeT(waypoints3',total_time);

waypoints4 = tubePaths(4).path;
nwaypoints = length(waypoints4(:,2));
distance = 0;
for k=2:nwaypoints
    distance = distance + norm(waypoints4(k,:) - waypoints4(k-1,:));
end
uavSpeed = 10;
total_time = distance / uavSpeed;
ts4 = arrangeT(waypoints4',total_time);

ts = (ts1 + ts2 + ts3 + ts4) / 4;


n_order = 6;
v0 = [0 0 0];
a0 = [0 0 0];
v1 = [0 0 0];
a1 = [0 0 0];

path_tmp = [path(2:end-1,1), path(2:end-1,4)];
polys_x = minimum_snap_single_axis_simple(waypoints1(:,1)',ts,n_order,v0(1),...
    a0(1),v1(1),a1(1),3,path_tmp);
path_tmp = [path(2:end-1,2), path(2:end-1,4)];
polys_y = minimum_snap_single_axis_simple(waypoints1(:,2)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
path_tmp = [path(2:end-1,3), path(2:end-1,4)];
polys_z = minimum_snap_single_axis_simple(waypoints1(:,3)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
num = 100;
gen_curve1 = GeneratorCalculate(polys_x,polys_y,polys_z,ts,num);
hold on;
scatter3(waypoints1(:,1),waypoints1(:,2),waypoints1(:,3),'red','filled');
plot3(gen_curve1.origin(1,:),gen_curve1.origin(2,:),gen_curve1.origin(3,:),'-r',LineWidth=2);

path_tmp = [path(:,1), path(:,4)];
polys_x = minimum_snap_single_axis_simple(waypoints2(:,1)',ts,n_order,v0(1),...
    a0(1),v1(1),a1(1),3,path_tmp);
path_tmp = [path(:,2), path(:,4)];
polys_y = minimum_snap_single_axis_simple(waypoints2(:,2)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
path_tmp = [path(:,3), path(:,4)];
polys_z = minimum_snap_single_axis_simple(waypoints2(:,3)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
gen_curve2 = GeneratorCalculate(polys_x,polys_y,polys_z,ts,num);
hold on;
scatter3(waypoints2(:,1),waypoints2(:,2),waypoints2(:,3),'red','filled');
plot3(gen_curve2.origin(1,:),gen_curve2.origin(2,:),gen_curve2.origin(3,:),'-r',LineWidth=2);

path_tmp = [path(:,1), path(:,4)];
polys_x = minimum_snap_single_axis_simple(waypoints3(:,1)',ts,n_order,v0(1),...
    a0(1),v1(1),a1(1),3,path_tmp);
path_tmp = [path(:,2), path(:,4)];
polys_y = minimum_snap_single_axis_simple(waypoints3(:,2)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
path_tmp = [path(:,3), path(:,4)];
polys_z = minimum_snap_single_axis_simple(waypoints3(:,3)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
gen_curve3 = GeneratorCalculate(polys_x,polys_y,polys_z,ts,num);
hold on;
scatter3(waypoints3(:,1),waypoints3(:,2),waypoints3(:,3),'red','filled');
plot3(gen_curve3.origin(1,:),gen_curve3.origin(2,:),gen_curve3.origin(3,:),'-r',LineWidth=2);

path_tmp = [path(:,1), path(:,4)];
polys_x = minimum_snap_single_axis_simple(waypoints4(:,1)',ts,n_order,v0(1),...
    a0(1),v1(1),a1(1),3,path_tmp);
path_tmp = [path(:,2), path(:,4)];
polys_y = minimum_snap_single_axis_simple(waypoints4(:,2)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
path_tmp = [path(:,3), path(:,4)];
polys_z = minimum_snap_single_axis_simple(waypoints4(:,3)',ts,n_order,v0(2),...
    a0(2),v1(2),a1(2),3,path_tmp);
gen_curve4 = GeneratorCalculate(polys_x,polys_y,polys_z,ts,num);
hold on;
scatter3(waypoints4(:,1),waypoints4(:,2),waypoints4(:,3),'red','filled');
plot3(gen_curve4.origin(1,:),gen_curve4.origin(2,:),gen_curve4.origin(3,:),'-r',LineWidth=2);

% 
% 


