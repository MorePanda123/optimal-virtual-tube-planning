# Optimal virtual tube planning and control for swarm robotics
This is the code for the paper: Mao, Pengda, Rao Fu, and Quan Quan. "Optimal virtual tube planning and control for swarm robotics." The International Journal of Robotics Research 43.5 (2024): 602-627. <https://doi.org/10.1177/02783649231210012>

# Introduce the framework of the codes
The codes are implemented in MATLAB, which are divided into two parts: 
1. Demo of optimal virtual tube planning in the random obstacled environments
2. Demo of MPC within the optimal virtual tube
These two parts code could be run independently and are without any dependances.

# How to run the codes
First, you should git clone the codes:
```
git clone https://github.com/MorePanda123/optimal-virtual-tube-planning.git
```
Then, enter the folder "optimal-virtual-tube-planning". Afer that, you will see two folders named "planning_3d" and "mpc" respectively. You can select any of them to run independtly, accroding to your usage.

## Demo of optimal virtual tube planning
Enter the folder "planning_3d" and run the file named "test.m" in MATLAB. Then, you will see the planning results in the figure.
## Demo of MPC within the optimal virtual tube
Enter the folder "mpc" and run the file named "main.m" in MATLAB. Then, you will see the results in the figure.

# Acknowledgement
The code of this paper refers to MPCC repo and Tube RRT* repo.