function [A,B,C,l1,l2] = state_space_inverted_pendulum()

    g = 9.81;     
    m1 = 0.5;      % mass of Arm 1 
    m2 = 0.3;      % mass of Arm 2 
    l1 = 0.3;      % length of Arm 1 
    l2 = 0.25;     % length of Arm 2 
    Ieq = 0.05;     % equivalent moment of inertia of pillar 
    
    %The above values chosen are completely arbitrary and serve only as a
    %conceptual bridge to visualise the state space equations.
    
    n = 0;
    phi_eq = n*pi;
    
    %This is the condition for stable/unstable equilibria.
    %From the derivation in the theory part:-
    %   phi = 0, 2*pi, 4*pi, ... -> Stable equilibria
    %   phi = pi, 3*pi, 5*pi, ... -> Unstable equilibria
    
    
    %Part 1: State Space Matrices
    
    D = l2*Ieq + (1/3)*m1*l1^2*l2 + (1/4)*m2*l1^2*l2;
    
    alpha_2_3 = (3/4)*(m2*g*l1*l2)/D;
    alpha_4_3 = (((-1)^(n+1)/2)*(3*g*Ieq + (m1+3*m2)*g*l1^2))/D;
    
    beta_2_1 = l2/D;
    beta_4_1 = (3*(-1)^(n+1)/2)*(l1/D);
    beta_2_2 = (3*(-1)^(n+1)/2)*(l1/D);
    beta_4_2 = (3*Ieq + (m1+3*m2)*l1^2)/(m2*l2*D);
    
    A = [0,  1,          0,  0;
         0,  0,  alpha_2_3,  0;
         0,  0,          0,  1;
         0,  0,  alpha_4_3,  0]; %State MAtrix
    
    B = [0,        0;
         beta_2_1, beta_2_2;
         0,        0;
         beta_4_1, beta_4_2]; %Control Matrix
    
    C = [1, 0, 0, 0;
         0, 0, 1, 0]; %Output Matrix

eig_numerical  = eig(A);
eig_analytical = [0; 0; sqrt(alpha_4_3); -sqrt(alpha_4_3)];

% The above has been directly implemented from the theory derivation. Next,
% we display all the values of parameters computed.

%
% 
% fprintf('D = %.6f\n', D);
% fprintf('alpha_2_3 = %.6f\n', alpha_2_3);
% fprintf('alpha_4_3 = %.6f\n', alpha_4_3);
% fprintf('beta_2_1 = %.6f\n', beta_2_1);
% fprintf('beta_4_1 = %.6f\n', beta_4_1);
% fprintf('beta_2_2 = %.6f\n', beta_2_2);
% fprintf('beta_4_2 = %.6f\n\n', beta_4_2);
% fprintf('A =\n'); 
% disp(A);
% fprintf('B =\n'); 
% disp(B);
% fprintf('C =\n'); 
% disp(C);
% fprintf('Numerical eigenvalues:\n'); disp(eig_numerical);
% fprintf('Analytical eigenvalues (Page 18):\n'); disp(eig_analytical);
% 
% if all(real(eig_numerical) <= 1e-10)
%     fprintf('=> Equilibrium phi = %d*pi is STABLE (n = %d)\n\n', n, n);
% else
%     fprintf('=> Equilibrium phi = %d*pi is UNSTABLE (n = %d)\n\n', n, n);
% end
% 
% 
% %Part 2: Transfer Function
% 
% omega = logspace(-2,3,1000);
% n_0 = length(omega);
% H = zeros(2,2,n_0);
% 
% for k = 1:n_0
%     s = 1j*omega(k);
%     H(:, :, k) = C*(inv(s*eye(4) - A))*B; %Required Transfer Function
% end
% 
% H11_mag = 20*log10(abs(squeeze(H(1,1,:))));
% H21_mag = 20*log10(abs(squeeze(H(2,1,:))));
% H12_mag = 20*log10(abs(squeeze(H(1,2,:))));
% H22_mag = 20*log10(abs(squeeze(H(2,2,:))));
% 
% H11_ph = angle(squeeze(H(1,1,:))) * 180/pi;
% H21_ph = angle(squeeze(H(2,1,:))) * 180/pi;
% H12_ph = angle(squeeze(H(1,2,:))) * 180/pi;
% H22_ph = angle(squeeze(H(2,2,:))) * 180/pi;
% 
% % Keeping the magnitude and phase separate to make plotting easier
% 
% 
% %Part 3: Plots
% 
% figure();
% 
% % a) Phase Portrait: We find the matrix exponential and use it to
% %    find the values of the state variables.
% 
% t = linspace(0, 10, 10000);
% x = zeros(length(t), 4);
% delta_x0 = [0.05; 0; 0.05; 0]; %Small perturbation from stable equilibrium
% 
% for k = 1:length(t)
%     x(k, :) = (expm(A*t(k)) * delta_x0)'; 
% end
% 
% figure(1);
% subplot(2,1,1);
% plot(x(:,1), x(:,2));
% xlabel('$\theta$ (rad)','Interpreter','latex')
% ylabel('$\dot{\theta}$ (rad/s)','Interpreter','latex')
% title('Phase Portrait 1');
% 
% subplot(2,1,2);
% plot(x(:,3), x(:,4));
% xlabel('$\phi$ (rad)','Interpreter','latex')
% ylabel('$\dot{\phi}$ (rad/s)','Interpreter','latex')
% title('Phase Portrait 2');
% 
% % b) Transfer Function Magnitude Response (|H(jw)|)
% 
% figure(2);
% 
% subplot(2,2,1);
% plot(log10(omega), H11_mag);
% xlabel('Frequency (rad/s)'); 
% ylabel('Magnitude (dB)');
% title('H_{11} magnitude');
% 
% subplot(2,2,2);
% plot(log10(omega), H12_mag);
% xlabel('Frequency (rad/s)'); 
% ylabel('Magnitude (dB)');
% title('H_{12} magnitude');
% 
% subplot(2,2,3);
% plot(log10(omega), H21_mag);
% xlabel('Frequency (rad/s)'); 
% ylabel('Magnitude (dB)');
% title('H_{21} magnitude');
% 
% subplot(2,2,4);
% plot(log10(omega), H22_mag);
% xlabel('Frequency (rad/s)'); 
% ylabel('Magnitude (dB)');
% title('H_{22} magnitude');
% 
% sgtitle('Transfer Function Magnitude Response');
% 
% % c) Transfer Function Phase Response
% 
% figure(3);
% 
% subplot(2,2,1);
% plot(log10(omega), H11_ph);
% xlabel('Frequency (rad/s)'); 
% ylabel('Magnitude (dB)');
% title('H_{11} phase');
% 
% subplot(2,2,2);
% plot(log10(omega), H12_ph);
% xlabel('Frequency (rad/s)'); 
% ylabel('Phase (deg)');
% title('H_{12} phase');
% 
% subplot(2,2,3);
% plot(log10(omega), H21_ph);
% xlabel('Frequency (rad/s)'); 
% ylabel('Phase (deg)');
% title('H_{21} phase');
% 
% subplot(2,2,4);
% plot(log10(omega), H22_ph);
% xlabel('Frequency (rad/s)'); 
% ylabel('Phase (deg)');
% title('H_{22} phase');
% 
% sgtitle('Transfer Function Phase Response');
% 
% disp('Script Completed Successfully!');

end