# Gyro-tower Simulation

> Modeling and Simulation in Science, Engineering, and Economics | Courant Institute of Mathematical Sciences | New York University

## Introduction

The gyroscope is simulated by the wheel constructed from networks of springs. More specifically, in each network, two points are connected with a dashpot and a spring, where the dashpots have their dashpot coefficients (the larger, the harder to pull) and the springs have their rest lengths and stiffness. The gyroscope is then constructed from several networks, attached to each other as in a bicycle wheel.

We know that a gyroscope is a mechanic system with a spinning wheel that harnesses the principle of conservation of angular momentum. This means that the spin of a the gyroscope tends to remain constant unless it is subjected to some external torque.

Then a natrual question may arise: what will happen when several gyroscopes are stacked up to form a **gyroscope tower**?

Well, an intial guess is that, if two gyroscopes are stacked together, each with opposite spinning direction, then such a stack may cancel out the gyroscopic effects. As long as the system holds together, if the two halves spin with exactly the same magnitude but opposite sign angular momentum, the system will behave just like a system of zero angular momentum.

However, from the standpoint of each rotating component, each requires a torque to change its own angular momentum. The angular momentums of the two separate components are to change extremely fast. The two spinning components must therefore exert huge torques on each other to achieve this.

In order to figure out what really happens, experiments will be done by modeling and statistical analysis of the gyroscope tower. Modeling and simulation will be done using MATLAB.

## Requirements

You need to have MATLAB installed. The programs are tested on version 9.12.0, 1855793. Other versions of MATLAB are not guaranteed to work.

1. Clone the GitHub repository to your local device:
```
git clone https://github.com/Charlie-XIAO/Gyro-simulation.git
```

2. Run `gyro_tower.m` in MATLAB, or use the following command line in the corresponding directory:
```
matlab -batch "gyro_tower; exit"
```

1. Note that when running `gyro_tower.m`, you will be prompted to input:

   - The number of gyroscopes in the gyroscope tower, positive integer expected.

   - Whether or not to make adjacent gyroscopes spin in opposite directions, `true` or `false` expected.

   - Whether or not to save a video of the simulation, `true` or `false` expected. If you are doing so, make sure you have `./videos` directory created.

You may figure out the code structure as follows:

```
Gyro-simulation
 │  gyro_tower.m        Main code for gyroscope tower simulation.
 |  wheel.m             The wheel function, returns data for modeling.
```

## Results

Some simulation videos are shown as follows:

<img src="./examples/2-rev-1.0-30.0.gif" width=45%> <img src="./examples/3-rev-1.0-30.0.gif" width=45%> <img src="./examples/4-rev-1.0-30.0.gif" width=45%> <img src="./examples/5-rev-1.0-30.0.gif" width=45%> <img src="./examples/6-rev-1.0-30.0.gif" width=45%> <img src="./examples/7-rev-1.0-30.0.gif" width=45%>

We have also done parameter analysis on the radius of the gyroscopes and their initial angular velocity. It turns out that the larger the radius and the larger the angular velocity, the more stable the gyroscope will become.

For more of our simulation and analysis results, check in [this folder](https://drive.google.com/drive/folders/1UpQE-VyRK9DaUnWq8wOW_s7JVZEpMOUg?usp=sharing) for our report, more simulation videos, and parameter analysis plots.

## References

Thanks for Professor [*Charles S. Peskin*](https://www.math.nyu.edu/~peskin/)'s tutoring and [notes](https://www.math.nyu.edu/~peskin/modsim_lecture_notes/index.html). Other references include [MATLAB](https://www.mathworks.com/help/matlab/) documentation, [Wikipedia page](https://en.wikipedia.org/wiki/Gyroscope) of gyroscopes, and Python library [pygyro](https://github.com/pyccel/pygyro). Also thanks Teaching Assistant [Mengjian Hua](mh5113@nyu.edu) for his suggestion on this project.
