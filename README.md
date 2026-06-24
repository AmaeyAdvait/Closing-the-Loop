# Closing the Loop
A Unified Approach to State-Space Kalman Decomposition, Reduced Order Luenberger Observer Design and Optimal LQR Control with the aim of modelling the Furata Pendulum

This is a project implemented as a part of the course: 'Systems and Control(EC31203)', IIT Kharagpur.

We start with an examination of controllability and observability using the Popov–Belevitch–Hautus (PBH) tests, followed by Kalman decomposition to isolate the intersections of controllable, uncontrollable and observable, unobservable subspaces of the system and use them to construct the similarity
transformation.

A reduced-order Luenberger observer is then constructed to estimate detectable states, ensuring asymptotic convergence of the estimation error. Next, an optimal state-feedback controller is designed using the Linear Quadratic Regulator (LQR) framework, balancing performance and control effort through the solution of the Discrete Algebraic Riccati Equation(DARE). These components are combined to form a complete output-feedback control law.

To validate the theoretical design, the closed-loop dynamics are analyzed through eigenvalue inspection and phase portrait visualization. The frame
work is further applied to a nonlinear benchmark system, the Furuta pendulum, where the effectiveness of the linearized control strategy is illustrated through animation. Finally, as an additional bonus on top of the mandatory requirement, the impact of stochastic disturbances is introduced through additive white Gaussian noise, providing insight into robustness and practical performance limitations.
