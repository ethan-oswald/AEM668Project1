%% project1_solution.m

%% Setup workspace.

close all;
clear all;
clc;


%% Set conversion constants.

RAD2DEG = 180/pi;
DEG2RAD = pi/180;


%% Set trim, geometric, mass, and aerodynamic parameters.
% All units are in meters, kilograms, Newtons, radians, and seconds.

% Calculate trim parameters.
trim.altitude = 6000;
trim.density = 0.66011;
trim.gravity = 9.788;
trim.speed_sound = 316.43;
trim.airspeed = 120;
trim.dynamic_pressure = 0.5*trim.density*trim.airspeed^2;
trim.mach = trim.airspeed/trim.speed_sound;
trim.flight_path = 0;
trim.coef_lift = 0.2984;

% Calculate the airplane parameters for the full vehicle.
full.coef_lift_0 = 0.01;
full.coef_drag_0 = 0.036;
full.coef_m_0 = 0.04;
full.coef_lift_slope = 5.05;
full.mass = 5800;
full.inertia_x = 112000;
full.inertia_y = 93000;
full.inertia_z = 194000;
full.gravity_center = 7.5;
full.aero_center = 6.8;

% Calculate wing parameters.
wing.area = 40;
wing.span = 20;
wing.tip_chord = 1;
wing.root_chord = 3;
wing.mean_chord = 2;
wing.x = full.gravity_center - full.aero_center;
wing.coef_lift_slope = 4.95;
wing.eff = 0.8;
wing.dihedral = 0.04;
wing.sweep = 0.1;
wing.coef_l_dihedral = -0.7;
wing.aspect_ratio = wing.span^2/wing.area;
wing.chord_ratio = wing.tip_chord/wing.root_chord;
wing.ail_inner = 7.5;
wing.ail_outer = 9.5;
wing.ail_factor = -0.1;
wing.ail_area = 0.35;
wing.ail_tau = interp1(0:0.05:0.7, [0, 0.16, 0.26, 0.34, 0.41, ...
    0.47, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.75, 0.78, 0.8], ...
    wing.ail_area/wing.area);

% Calculate vertical tail parameters.
vtail.area = 5;
vtail.span = 2.5;
vtail.mean_chord = 2;
vtail.x = -8.5;
vtail.z = -0.8;
vtail.coef_lift_slope = 3;
vtail.eff = 0.95;
vtail.sidewash_slope = 0.1;
vtail.volume_ratio = -vtail.x*vtail.area/(wing.area*wing.span);
vtail.rudder_area = 1.2;
vtail.rudder_tau = interp1(0:0.05:0.7, [0, 0.16, 0.26, 0.34, 0.41, ...
    0.47, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.75, 0.78, 0.8], ...
    vtail.rudder_area/vtail.area);


%Derivation of Stability and Control Derivatives
der.C_Y_beta = -vtail.eff * (vtail.area/wing.area) * vtail.coef_lift_slope * (1 + vtail.sidewash_slope);
der.Y_beta = (trim.dynamic_pressure * wing.area)/full.mass * der.C_Y_beta;
der.C_L_beta = wing.coef_l_dihedral * wing.dihedral;
der.L_beta = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_x * der.C_L_beta;
der.C_N_beta = 0 + 0.95 * vtail.volume_ratio * vtail.coef_lift_slope * (1 + vtail.sidewash_slope);
der.N_beta = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_z * der.C_N_beta;
der.C_Y_p = (wing.aspect_ratio + cos(wing.sweep))/(wing.aspect_ratio + 4 * cos(wing.sweep)) * tan(wing.sweep) * trim.coef_lift;
der.Y_p = (trim.dynamic_pressure * wing.area * wing.span)/(2 * full.mass * trim.airspeed) * der.C_Y_p; 
der.C_L_p = -(1 + 3*wing.chord_ratio)/(1+wing.chord_ratio) * wing.coef_lift_slope/12;
der.L_p = (trim.dynamic_pressure * wing.area * wing.span^2)/(2*full.inertia_x*trim.airspeed)*der.C_L_p;
der.C_N_p = -trim.coef_lift/8;
der.N_p = (trim.dynamic_pressure * wing.area * wing.span^2)/(2*full.inertia_z*trim.airspeed)*der.C_N_p;
der.C_Y_r = 2 * vtail.x/wing.span * der.C_Y_beta;
der.Y_r = (trim.dynamic_pressure * wing.area * wing.span)/(2*full.mass*trim.airspeed) * der.C_Y_r;
der.C_L_r = trim.coef_lift/4 - 2*(vtail.x*vtail.z/wing.span^2)*der.C_Y_beta;
der.L_r = (trim.dynamic_pressure * wing.area * wing.span^2)/(2*full.inertia_x*trim.airspeed) * der.C_L_r;
der.C_N_r = 2 * vtail.eff * vtail.volume_ratio * vtail.x/wing.span * vtail.coef_lift_slope;
der.N_r = (trim.dynamic_pressure * wing.area * wing.span^2)/(2*full.inertia_z * trim.airspeed) * der.C_N_r;
der.C_Y_delta_r = vtail.area/wing.area * vtail.rudder_tau * vtail.coef_lift_slope;
der.Y_delta_r = (trim.dynamic_pressure * wing.area)/full.mass * der.C_Y_delta_r;
der.C_L_delta_r = vtail.area/wing.area * vtail.rudder_tau * vtail.z/wing.span * vtail.coef_lift_slope;
der.L_delta_r = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_x * der.C_L_delta_r;
der.C_N_delta_r = -vtail.eff * vtail.volume_ratio * vtail.coef_lift_slope * vtail.rudder_tau;
der.N_delta_r = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_z * der.C_N_delta_r;
der.C_L_delta_a = (2 * wing.coef_lift_slope * wing.ail_tau * wing.root_chord)/(wing.area * wing.span) * ((wing.ail_outer^2/2 + (wing.tip_chord/wing.root_chord - 1)/(wing.span/2) * wing.ail_outer^3/3)-(wing.ail_inner^2/2 + (wing.tip_chord/wing.root_chord - 1)/(wing.span/2) * wing.ail_inner^3/3));
der.L_delta_a = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_x * der.C_L_delta_a;
der.C_N_delta_a = 2 * wing.ail_factor * trim.coef_lift * der.C_L_delta_a;
der.N_delta_a = (trim.dynamic_pressure * wing.area * wing.span)/full.inertia_z * der.C_N_delta_a;

K_ari = -der.N_delta_a/der.N_delta_r;


%% Calculation of Lateral-Directional Matrices
A_lat = zeros(4,4);
A_lat(1,1) = der.Y_beta/trim.airspeed;
A_lat(1,2) = der.Y_p/trim.airspeed;
A_lat(1,3) = der.Y_r/trim.airspeed - 1;
A_lat(1,4) = trim.gravity/trim.airspeed;
A_lat(2,1) = der.L_beta;
A_lat(2,2) = der.L_p;
A_lat(2,3) = der.L_r;
A_lat(3,1) = der.N_beta;
A_lat(3,2) = der.N_p;
A_lat(3,3) = der.N_r;
A_lat(4,2) = 1;

B_lat = zeros(4,2);
B_lat(1,2) = der.Y_delta_r / trim.airspeed;
B_lat(2,1) = der.L_delta_a;
B_lat(2,2) = der.L_delta_r;
B_lat(3,1) = der.N_delta_a;
B_lat(3,2) = der.N_delta_r;

%Print Answer for Part a.

fprintf("The 4x4 lateral-directional state matrix is\n");
display(A_lat);
fprintf("The 4x2 input matrix is \n")
display(B_lat);




%% Calculation of Lateral-Directional Modes
eigenvalues = eig(A_lat);
%Dutch Roll Mode
DR.quad = [1, (-der.N_r - der.Y_beta/trim.airspeed), (der.N_beta + der.N_r * der.Y_beta/trim.airspeed)];
DR.eig = eigenvalues(2);
DR.freq = sqrt(DR.quad(3));
DR.damp = (DR.quad(2)/(2*DR.freq));

%Roll Mode
roll.quad = [1, -der.L_p];
roll.eig = eigenvalues(3);


%Spiral Mode 
spiral.quad = [1, trim.gravity * (der.L_beta * der.N_r - der.L_r * der.N_beta)/(der.Y_beta * (der.L_r*der.N_p - der.N_r*der.L_p) + trim.airspeed*(der.L_beta*der.N_p - der.N_beta*der.L_p))];
spiral.eig = eigenvalues(4);

%Print Answer for Part B
fprintf("The eigenvalues for the dutch roll mode are:\n")
display(eigenvalues(1))
display(eigenvalues(2))
fprintf("The natural frequency for the dutch roll mode is %5.2f\n", DR.freq);
fprintf("The damping ratio for the dutch roll mode is %5.2f\n", DR.damp);

fprintf("The eigenvalue for the roll mode is")
display(eigenvalues(3));
fprintf("The eigenvalue for the spiral mode is")
display(eigenvalues(4));
fprintf("Since the spiral mode is not in the LHP, the aircraft is not laterally-directionally stable.\n" + ...
    "Therefore, a yaw damper must be included to ensure stability.\n")


%% Creation of outer loop guidance system
plant = tf([0 trim.gravity], [trim.airspeed 0]);
loop_bandwidth = 0.8;

%Lead stage
freq_bar = loop_bandwidth;
beta = 2.2;
lead_stage = tf([beta freq_bar], [1 beta*freq_bar]);

%Integral Boost Stage
freq_bar = 2.4;
int_boost = tf([1 freq_bar], [1 0]);

%Rolloff stage
freq_bar = 5;
rolloff_stage = tf(freq_bar, [1 freq_bar]);

K = int_boost * lead_stage * rolloff_stage;

%Proportional Gain stage
loop_bandwidth = 0.8;
prop_stage = 1/abs(freqresp(plant*K,loop_bandwidth));

%Calculate and display requested values
outer.controller = K * prop_stage;
outer.open_loop = plant * K * prop_stage;
outer.closed_loop = feedback(outer.open_loop,1);
outer.poles = pole(outer.closed_loop);
[outer.gm, outer.pm] = margin(outer.open_loop);
outer.crossover = getGainCrossover(outer.open_loop,1);
outer.sensitivity = 1/(1+outer.open_loop);
outer.steady_state_error = freqresp(outer.sensitivity,0);
outer.sensitivity = frd(feedback(1, outer.open_loop), 0.05);
outer.comp_sensitivity = frd(feedback(outer.open_loop,1), 10);


fprintf("The heading hold transfer function is\n")
display(outer.closed_loop);
fprintf("The closed loop poles for the heading hold outer-loop are:\n ")
display(outer.poles);
fprintf("so stability conditions are met.\n")
fprintf("The gain margin for the outer-loop is %5.2f dB\n", 20*log10(outer.gm));
fprintf("The phase margin for the outer-loop is %5.2f degrees\n", outer.pm);
fprintf("The loop bandwidth is %5.2f rad/s\n", outer.crossover);
fprintf("The steady state tracking error for unit step inputs is %5.2f\n", outer.steady_state_error);
fprintf("The tracking error for frequencies below 0.05 rad/s is <= %4.2f percent\n", 100*abs(outer.sensitivity.Response));
fprintf("The gain at frequencies above 10 rad/s is <= %4.2f percent\n\n\n", 100*abs(outer.comp_sensitivity.Response))

%% Inner Loop

%Creation of C and D Matrices for roll controllers
C = [0 0 0 1; 0 0 1 0];
D = [0 0;0 0];

%Creation of Actuator tf and Matrix
actuator = tf(100, [1, 100]);
Act_mat = [actuator, 0; 0, actuator];

%Creation of Yaw Damper
K_r = tf([1 0], [4.2 .01]);

G_ss = ss(A_lat, B_lat, C, D);
G = tf(G_ss) * Act_mat;
inner.tf1 = G(1,1);

inner.loop_bandwidth = 10;

%Creation of controller for first transfer function
%Integral Boost stage
freq_bar = 5;
integral_boost = tf([1 freq_bar], [1 0]);
%Lead stage
beta = 5;
lead_stage = tf([beta inner.loop_bandwidth], [1 beta*inner.loop_bandwidth]);
K = integral_boost * lead_stage;
%proportional gain stage
prop_stage = 1/abs(freqresp(inner.tf1*K,inner.loop_bandwidth));


inner.controller = K * prop_stage;
inner.openloop1 = inner.tf1 * inner.controller;

%Calculate and display requested values
inner.closed_loop1 = feedback(inner.openloop1,1);
inner.poles1 = pole(inner.closed_loop1);
[inner.gm1, inner.pm1] = margin(inner.openloop1);
inner.crossover1 = getGainCrossover(inner.openloop1,1);
inner.sensitivity1 = 1/(1+inner.openloop1);
inner.steady_state_error1 = freqresp(inner.sensitivity1,0);
inner.sensitivity1 = frd(feedback(1, inner.openloop1), 1);
inner.error1 = 100*abs(inner.sensitivity1.Response);
inner.comp_sensitivity1 = frd(feedback(inner.openloop1,1), 100);
inner.gain1 = 100 * abs(inner.comp_sensitivity1.response);

fprintf("The roll control transfer function is\n")
display(inner.closed_loop1);
fprintf("The transfer function for the rudder and aileron actuators are\n")
display(actuator);
fprintf("The transfer function for the yaw damper is\n")
display(K_r);
fprintf("The closed loop poles for the heading hold outer-loop are:\n ")
display(inner.poles1);
fprintf("so stability conditions are met.\n")
fprintf("The gain margin for the inner-loop is %5.2f dB\n", 20*log10(inner.gm1));
fprintf("The phase margin for the inner-loop is %5.2f degrees\n", inner.pm1);
fprintf("The loop bandwidth is %5.2f rad/s\n", inner.crossover1);
fprintf("The steady state tracking error for unit step inputs is %5.2f\n", inner.steady_state_error1);
fprintf("The tracking error for frequencies below 1 rad/s is <= %4.2f percent\n", 100*abs(inner.sensitivity1.Response));
fprintf("The gain at frequencies above 100 rad/s is <= %4.2f percent\n", 100*abs(inner.comp_sensitivity1.Response));


%% Run Simulink simulation and retrieve values
sim = sim("OswaldEProj1Sim");
open('OswaldEProj1Sim')
headingseries = sim.heading;
time = headingseries.time;
heading = headingseries.data;

%Plot heading and output requested values
plot(headingseries);
ylabel("Heading angle(degrees)")
S = stepinfo(heading,time);
fprintf("\nThe simulated settling time is %5.2f seconds\n", S.SettlingTime);
fprintf("\nThe Maximum Overshoot is %5.2f percent\n", S.Overshoot);



