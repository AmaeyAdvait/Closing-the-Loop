clc;
clear;
close all;

figure;


% Programmer's Log: The current and complete implementation takes 
% approximately 3 minutes to run on a standard personal computer. 
% This file is the main file where the simulation can be generated. The 
% file generates 5 figures and one animation that runs for a duration of 
% approximately 1 minute. This file spans 441 lines of code. Apart from 
% this main file, the total project consists of 5 helper function files:-
%
% 1. gram_schmidt.m : Computes the orthonormal basis vectors of a given 
% matrix using gram schmidt orthogonalization. It spans 16 lines of code.
%
% 2. Kalman_Decomposition.m: Computes the similarity transformation T using
% the Zassenhaus Intersection method. It spans 151 lines of code.
%
% 3. Luenberger_Observer.m: Computes the observer gain matrix L using the
% observability Gramian matrix and verifying positive semi-definiteness for
% Lyapunov Stability. It spans 71 lines of code.
%
% 4. LQR.m: Computes the state gain matrix K by solving the Algebraic 
% Ricatti Equation using the Hamiltonian matrix and directly applying the
% obtained value to u =-kx. It spans 23 lines of code.
%
% 5. state_space_inverted_pendulum.m: Computes the required A,B,C (state
% matrix, control input matrix and output matrix) for the double inverted 
% pendulum and passes the lengths of the pendulum rods(arbitrarily chosen). 
% Directly obtained from assignment-1 with no changes required. 
% It spans 192 lines of code.
%
% A bonus section verifying robustness has also been implemented. A
% detailed guide regarding its formulation and implementation can be found
% on line 150.


function x_co = closing_the_loop(A,B,C,figure_count)

    %Part 1: Contollability and Observability Tests

    Y = eig(A);
    % disp(Y);
    N_A = length(Y);
    I_N_A = eye(N_A);
    
    %Controllability
    param_C = 0;
    for i = 1:length(Y)
        lambda = Y(i);
        PBH_matrix_C = [(lambda*I_N_A - A), B];
    
        if rank(PBH_matrix_C) < N_A
            param_C = 1;
        end
    end
    
    if param_C == 0
        disp('(A,B) is controllable');
    else
        disp('(A,B) is uncontrollable');
    end
    
    %Observability
    param_O = 0;
    for i = 1:length(Y)
        lambda = Y(i);
        PBH_matrix_O = [(A - lambda*I_N_A); C];
    
        if rank(PBH_matrix_O) < N_A
            param_O = 1;
        end
    end
    
    if param_O == 0
        disp('(C,A) is observable');
    else
        disp('(C,A) is unobservable');
    end
    

    % Part-2: Kalman Decomposition 
    Cm = B;
    for i = 1:N_A-1
        Cm = [Cm, A^i*B];
    end
    
    Om = C;
    for i = 1:N_A-1
        Om = [Om; C*A^i];
    end
    
    n_Cm = rank(Cm);
    n_Om = rank(Om);
    
    disp(n_Cm);
    disp(n_Om);
    
    [T, n_v_c_o, n_v_c_obar, n_v_cbar_o, n_v_cbar_obar] = Kalman_Decomposition(Cm,Om);
    disp(T);
    
    A_bar = T\(A*T);
    B_bar = T\B;
    C_bar = C*T;
    
    disp(A_bar);
    disp(B_bar);
    disp(C_bar);
    

    %Part-3: Reduced Order Luenberger Observer
    %Matrix formation: Picking out the observable blocks
    
    A11 = A_bar(n_v_c_obar+1:n_v_c_obar+n_v_c_o, n_v_c_obar+1:n_v_c_obar+n_v_c_o);
    A12 = A_bar(n_v_c_obar+1:n_v_c_obar+n_v_c_o, ...
        n_v_c_obar+n_v_c_o+n_v_cbar_obar+1:n_v_c_obar+n_v_c_o+n_v_cbar_obar+n_v_cbar_o);
    A21 = A_bar(n_v_c_obar+n_v_c_o+n_v_cbar_obar+1:n_v_c_obar+n_v_c_o+n_v_cbar_obar + n_v_cbar_o, ...
        n_v_c_obar+1:n_v_c_obar+n_v_c_o);
    A22 = A_bar(n_v_c_obar+n_v_c_o+n_v_cbar_obar+1:n_v_c_obar+n_v_c_o+n_v_cbar_obar+n_v_cbar_o, ...
        n_v_c_obar+n_v_c_o+n_v_cbar_obar+1:n_v_c_obar+n_v_c_o+n_v_cbar_obar+n_v_cbar_o);
    
    B1 = B_bar(n_v_c_obar+1:n_v_c_obar+n_v_c_o,:);
    B2 = B_bar(n_v_c_obar+n_v_c_o+n_v_cbar_obar+1:n_v_c_obar+n_v_c_o+n_v_cbar_obar+n_v_cbar_o, :);
    
    disp(A12);
    disp(A22);
    L = Luenberger_Observer(A12, A22);
    

    %Part-4: Finding and solving the Algebraic Ricatti Equation
    A_co = A_bar(n_v_c_obar+1:n_v_c_obar+n_v_c_o, n_v_c_obar+1:n_v_c_obar+n_v_c_o);
    B_co = B_bar(n_v_c_obar+1:n_v_c_obar+n_v_c_o, :);
    C_co = C_bar(:, n_v_c_obar+1:n_v_c_obar+n_v_c_o);
    
    K = LQR(A_co, B_co, C_co, n_v_c_o);
    

    %Designing the Luenberger-LQR(LL) closed loop
    if figure_count == 2
        A_cl = A_co - B_co*K - L*C_co; %Combined LQR+ROLO in the closed loop
    else
        A_cl = A_co - B_co*K; %Only LQR block in the closed loop
    end
    

    % Bonus: Checking for Robustness in the system

    %Guide to running this block: The threshold for the weight of the noise
    %(sigma) = 0.1. For sigma <= 0.1, the LQR control worked reliably 
    %across multiple test runs. In order to reliably observe test cases 
    %where LQR fails to control the system, set sigma = 1.

    % Programmer's Log: We introduce an additive white Gaussian noise and 
    % make it state dependent(linearly proportional in the base case) to 
    % observe the robustness of the system to noise. This was done as a
    % logical continuation to the assignment as explained in class to model
    % practical applicability of the controller.

    sigma = 0.1;
    eta = randn(size(A_cl,2));
    A_cl = A_cl + sigma*eta; % This has been modeled as a state-dependent multiplicative
                             % noise, more details in the report

    disp(' ');
    disp('Eigenvalues of the closed loop state matrix(A_cl)')
    disp(eig(A_cl));
    

    %Part-5: Plotting the closed-loop phase portrait
    t = linspace(0,4,1600);
    N_t = length(t);
    colour_count = 20;
    x_co = zeros(n_v_c_o, N_t);
    y_co = zeros(size(C_co,1), N_t);
    h = zeros(1,colour_count);
    colours = lines(colour_count);
    
    figure(figure_count);
    hold on;
    
    for k = 1:colour_count
        x_initial = randn(N_A,1);
        x_initial_t = T\x_initial; % Similarity transformation of initial input
        x_initial_co = x_initial_t(n_v_c_obar+1:n_v_c_obar+n_v_c_o,1);
    
        for i = 1:length(t)
            x_co(:,i) = expm((A_cl)*t(i))*x_initial_co;
            y_co(:,i,k) = real(C_co*x_co(:,i));
        end
        plot(y_co(1,:,k), y_co(2,:,k), 'Color', colours(k,:), 'lineWidth', 1.5);
        h(k) = plot(y_co(1,1,k), y_co(2,1,k), 'o', 'MarkerFaceColor', colours(k,:), ...
                   'MarkerSize', 6);
    end
    
    for i = 1:N_t
        for k = 1:colour_count
            set(h(k), 'XData', y_co(1,i,k), 'YData', y_co(2,i,k));
        end
        drawnow limitrate;
    end
    hold off;
    xline(0, 'k', 'LineWidth', 1);
    yline(0, 'k', 'LineWidth', 1);
    grid on;
    xlabel('y1');
    ylabel('y2');

    if figure_count == 2
        title('Output Phase Portrait of Problem-1: y1 vs y2');
    else
        title('Output Phase Portrait of Problem-2: y1 vs y2');
    end

end


%Problem-1: 

A = [-1.213213255856908  -0.052699540857194   0.114452537819400  -0.458503457153880  -0.147596519643947   0.330216171993562  -0.216306142694950   0.108524691675982  -0.277136364174417   0.213306124705031 -0.234400418544842  -0.694904339299891  -0.372985044024789  -0.121052395899524;
   -0.445138204194362  -1.087338372330932   0.506316210864373   0.420119557543333   0.489063639965910  -0.265671000243473  -0.336650153174219  -0.266447453528354   0.326773051102463  -0.796414105574600   -0.152486416125726  -0.360052673382867  -0.509891627577762   0.880961088466860;
   0.484525786600994   0.174803353899379  -0.108630387296726  -0.147473032367577  -0.370609448513368   0.123445036090085  -0.391828220999852   0.432403506843010   0.040561080322674   1.548617577592084   -0.375626137358579  -0.964927403783107   0.116685942887537   0.195907069263956;
  -0.628653100827204   0.441686129222833  -0.196853639475424  -0.740518371993041  -0.203043754424680  -0.312951739221684   0.459277316310521  -1.144480648269924  -0.227661184417581   0.621646540385829    0.773529276032394  -0.824173422452114  -0.511336138546041  -0.143886170459940;
  -0.526879261134897   0.496721354898758  -0.789140781845308  -0.351413813417401  -1.781412350680168   0.300097072779868   0.361009179150026   0.821998420305273   0.197670222666982   0.103894550459810   -0.518715568276280   0.458711846080412  -0.050433972884675   0.017439085643547;
  -0.162938036835619   0.115227089326697   0.230518490031468   0.226985493452003  -0.147401905420421  -1.073949490350060   0.492057044397738   0.809216490217785  -0.106022553236031  -0.200408629363208    0.219266594012242  -1.284874915967967  -0.219452186333267   0.086424084364871;
  -0.937922024260352   0.295340626262125   0.081848014134194   0.702150763890030   0.460918049982761  -0.410290741684858  -1.656970404290512   0.515521785489352   0.039634647589333   0.551677571844913    0.910747253955778  -0.502106754099405  -1.258898348872255   0.174292978828746;
   0.727293350401215  -0.399923820740790   0.252644270324252  -1.022331509304516   0.714779970671842   0.160297998212142  -0.683617227260968  -0.862980196841881  -0.434983508111977  -0.749659627595721   -0.848332491977190  -0.073701377614643   0.047794275307910  -0.291037548459371;
  -0.883478418688157   0.807423337770903  -0.272683699489586   0.267422703430696  -0.335917976087635  -0.355488407278101   1.117546803497625   0.475823305013847  -1.343243368591083   1.216242923480291   -0.117171815598981  -0.380364549164507  -0.827852126523950   0.367288610907304;
  -0.432084382768220  -0.208121861637220   1.692305102882171   0.996853318049892  -0.253373736345148  -1.076336539502446   0.349952751653385   0.243525258046507   0.467699882747476  -0.138850036162691    0.058652839684074  -0.660320581520629   0.238591279631722  -0.021125539243573;
  -0.162352919624227  -0.976155299865440  -0.243225326820792   0.836239633325489  -0.433802463517382   0.411365862850318   0.130321781352260  -0.897920092095636  -0.228424553357463  -0.568131183985648   -0.593732855151864  -0.725643142274004   0.182526229046144  -0.147385751504529;
  -0.264370906345628  -0.184742644441694  -1.003940450295694  -0.848369091986984   0.494546791639666  -1.104312423847451  -0.472200324956500  -0.339323185019347  -0.108101954544411  -0.430160178902030   -0.930981520180108   0.321055052198516   0.688882757485099  -0.611757590647088;
   0.575439200692222  -1.267310249593962  -0.215460248084041  -0.312938157224291   0.021629740368688   0.430071077232889  -1.163510493607768  -0.659695990919597  -0.382659856177330   0.430530215423247   -0.196984819655413   0.942614992066365  -0.194564541005406  -0.908760104887207;
  -0.214375145115562   0.987648714149219   0.442806720234774  -0.081107514151397  -0.216286104227307  -0.199881487887741   0.281210418831350  -0.068571314905706   0.190937766594154  -0.003347907463792   -0.101054426075363  -0.641224659163822  -1.242848391866568  -1.725651421647249];

B= [-0.051723763561734  -0.305270449184583  -0.015315868248594;
  -0.116732420024660  -0.475945089926195   0.149955878770381;
  -0.112898591360926  -0.041230352019466  -0.257329596734228;
   0.450703580365998  -0.195616995582868   0.021668047239566;
  -0.175175023018630  -0.164833619607650   0.126111050083874;
   0.193325591568345  -0.201825450079648   0.001812952706111;
  -0.746946115214500  -0.888561633091335  -0.534882101127978;
   0.118896143573814   0.003278035071589   0.438395454885899;
  -0.834413562808927  -0.496999682104461  -0.647097984062711;
  -0.486371252402602  -0.307207657585152  -0.536947426548773;
   0.293031461319545   0.647280729575534  -0.241442432395560;
   0.281522419773345   0.703605033443511   0.496061513668761;
   0.401731276643778   0.718304876069954   0.425248352781694;
  -0.207634639546252  -0.326879707471223  -0.260613735300743];

C =[0.695053806816924  -0.653901035211257  -0.216669043964540  -0.406643748943393  -0.337496928708196   0.800344919841996  -0.734680898447967  -0.383177337555783  -0.062081713564930  -0.164440009729601  -0.052596268307648   0.283621895206324   0.826159797514201  -0.281394936691068;
   0.170497399997020  -0.095533189106804  -0.211894215368707  -0.103936550394383  -0.119769092888598   0.354854849584245  -0.625481232333294  -0.485883175434064   0.198880592962716  -0.214578324733530    0.086360044518748   0.383477729131322   0.209925692065307   0.087002975029789];

Y = eig(A);
% disp(Y);
N_A = length(Y);

param_A = 0;
for i = 1:N_A
    if Y(i) > 0
        param_A = 1;
    end
end

disp('Problem-1');
disp(' ');
if param_A == 0
    disp('A is stable');
else
    disp('A is not stable');
end

t = linspace(0,0.5,1000);
N_t = length(t);
colour_count = 20;
x = zeros(N_A, N_t);
y = zeros(2, N_t);
h = zeros(1,colour_count);
colours = lines(colour_count);

figure(1);
hold on;

for k = 1:colour_count
    x_initial = randn(N_A,1);
    for i = 1:length(t)
        x(:,i) = expm(A*t(i))*x_initial;
        y(:,i,k) = C*x(:,i);
    end
    plot(y(1,:,k), y(2,:,k), 'Color', colours(k,:), 'lineWidth', 1.5);
    h(k) = plot(y(1,1,k), y(2,1,k), 'o', 'MarkerFaceColor', colours(k,:), 'MarkerSize', 6);
end

for i = 1:N_t
    for k = 1:colour_count
        set(h(k), 'XData', y(1,i,k), 'YData', y(2,i,k));
    end
    drawnow limitrate;
end
hold off;
xline(0, 'k', 'LineWidth', 1);
yline(0, 'k', 'LineWidth', 1);
grid on;
xlabel('y1');
ylabel('y2');
title('Input Phase Portrait of Problem-1: y1 vs y2');


figure_count = 2;
closing_the_loop(A,B,C,figure_count);
figure_count = figure_count+2;


%Problem-2: Designing an LQR for the Furuta Pendulum

[A,B,C,l1,l2] = state_space_inverted_pendulum();

Y = eig(A);
% disp(Y);
N_A = length(Y);

param_A = 0;
for i = 1:N_A
    if Y(i) > 0
        param_A = 1;
    end
end

disp(' ');
disp('Problem-2');
disp(' ');
if param_A == 0
    disp('A is stable');
else
    disp('A is not stable');
end

t = linspace(0,0.5,1000);
N_t = length(t);
colour_count = 20;
x = zeros(N_A, N_t);
y = zeros(2, N_t);
h = zeros(1,colour_count);
colours = lines(colour_count);

figure(3);
hold on;

for k = 1:colour_count
    x_initial = randn(N_A,1);
    for i = 1:length(t)
        x(:,i) = expm(A*t(i))*x_initial;
        y(:,i,k) = C*x(:,i);
    end
    plot(y(1,:,k), y(2,:,k), 'Color', colours(k,:), 'lineWidth', 1.5);
    h(k) = plot(y(1,1,k), y(2,1,k), 'o', 'MarkerFaceColor', colours(k,:), 'MarkerSize', 6);
end

for i = 1:N_t
    for k = 1:colour_count
        set(h(k), 'XData', y(1,i,k), 'YData', y(2,i,k));
    end
    drawnow limitrate;
end
hold off;
xline(0, 'k', 'LineWidth', 1);
yline(0, 'k', 'LineWidth', 1);
grid on;
xlabel('y1');
ylabel('y2');
title('Input Phase Portrait of Problem-2: y1 vs y2');

x_p = closing_the_loop(A,B,C,figure_count);
x_p = real(x_p);

%Animation
figure('Name','Furuta Pendulum Animation','NumberTitle','off', ...
    'Position',[100 100 700 500]);

t_p = linspace(0,100,1600);

for ti = 1:3:N_t

    theta = x_p(1,ti);   % base rotation
    phi   = pi + x_p(3,ti);   % pendulum angle

    % Arm
    x_base = l1*cos(theta);
    y_base = l1*sin(theta);
    z_base = 0;

    % Pendulum tip
    px = x_base + l2*sin(phi)*cos(theta);
    py = y_base + l2*sin(phi)*sin(theta);
    pz = -l2*cos(phi);

    clf;
    figure(gcf);

    subplot(1,2,1);
    hold on;

    % Rod 1 (base arm)
    plot3([0, x_base], [0, y_base], [0, z_base], ...
        'b','LineWidth',3);

    % Rod 2 (pendulum)
    plot3([x_base, px], [y_base, py], [z_base, pz], ...
        'k','LineWidth',3);

    axis equal;
    axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);
    grid on;
    view(45,25);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('t = %.2f s', t_p(ti)));

    subplot(1,2,2);
    plot(t_p(1:ti), x_p(3,1:ti)*(180/pi),'b','LineWidth',1.5); hold on;
    plot(t_p(1:ti), x_p(1,1:ti),'r','LineWidth',1.5);

    xlabel('Time (s)');
    legend('\phi (deg)', '\theta (rad)','Location','northeast');
    title('State trajectories');
    grid on;

    drawnow;
end

fprintf('Animation complete.\n');

figure(5);
subplot(2,2,1); plot(t_p, x_p(1,:)); grid on; xlabel('t'); ylabel('s (m)'); title('\theta');
subplot(2,2,2); plot(t_p, x_p(2,:)); grid on; xlabel('t'); ylabel('ds/dt'); title('$\dot{\theta}$', 'Interpreter', 'latex');
subplot(2,2,3); plot(t_p, x_p(3,:)*180/pi); grid on; xlabel('t'); ylabel('\theta (deg)'); title('\phi');
subplot(2,2,4); plot(t_p, x_p(4,:)); grid on; xlabel('t'); ylabel('d\theta/dt'); title('$\dot{\phi}$', 'Interpreter', 'latex');
sgtitle('Inverted Pendulum — LQR Closed Loop trajectories');

disp(' ');
disp('Completed Successfully!');