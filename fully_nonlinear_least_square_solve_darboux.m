function [u_k,v_k]=fully_nonlinear_least_square_solve_darboux(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight)
prob = optimproblem;
% options = optimoptions('fminunc','Display','none');
% options = optimoptions('fminunc','Display','iter-detailed','FunctionTolerance',1e-6,'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
options = optimoptions('fminunc','Display','none','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'StepTolerance',1e-12);
helperfun=@(u,v) nonlinear_least_square_solve_darboux_helper(R_helix,R_straight,omega_helix,omega_straight,theta_helix_0,theta_straight_0,n_helix,n_straight,s,u0_k,v0_k,epsilon_p_helix,epsilon_p_straight,u,v);
if n_helix+n_straight>=6
    u_k = optimvar('u_k',3,1);
    v_k = optimvar('v_k',3,1);
    prob.Objective=helperfun([u_k(1);u_k(2);u_k(3)],[v_k(1);v_k(2);v_k(3)]);
    sol0.u_k=u0_k;
    sol0.v_k=v0_k;
    % helperfun([0;0;0],[0;0;1])
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=sol.u_k;
    v_k=sol.v_k;
elseif n_helix+n_straight==5
    u_k = optimvar('u_k',2,1);
    v_k = optimvar('v_k',3,1);
    prob.Objective=helperfun([u_k(1);u_k(2);u0_k(3)],[v_k(1);v_k(2);v_k(3)]);
    sol0.u_k=[u0_k(1);u0_k(2)];
    sol0.v_k=v0_k;
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=[sol.u_k;u0_k(3)];
    v_k=sol.v_k;
elseif n_helix+n_straight==4
    % u_k = optimvar('u_k',3,1);
    % v_k = optimvar('v_k',1,1);
    % prob.Objective=helperfun([u_k(1);u_k(2);u_k(3)],[v0_k(1);v0_k(2);v_k(1)]);
    % sol0.u_k=u0_k;
    % sol0.v_k=[v0_k(3)];
    % prob.Constraints.cons1=u_k(3)<=0.00001;
    % prob.Constraints.cons2=u_k(3)>=-0.00001;
    % % prob.Constraints.cons3=u_k(2)<=10;
    % % prob.Constraints.cons4=u_k(2)>=-10;
    % % prob.Constraints.cons3=u_k(1)<=10;
    % % prob.Constraints.cons4=u_k(1)>=-10;
    % % prob.Constraints.cons3=v_k(1)<=0.8;
    % % prob.Constraints.cons4=v_k(1)>=1.2;
    % options = optimoptions('fmincon','Display','none','FunctionTolerance',1e-20,'OptimalityTolerance',1e-20,'StepTolerance',1e-20);
    % sol = solve(prob,sol0,'Options',options,'Solver','fmincon');
    % u_k=sol.u_k;
    % v_k=[v0_k(1);v0_k(2);sol.v_k];

    u_k = optimvar('u_k',2,1);
    v_k = optimvar('v_k',1,1);
    prob.Objective=helperfun([u_k(1);u_k(2);u0_k(3)],[v0_k(1);v0_k(2);v_k(1)]);
    sol0.u_k=[u0_k(1);u0_k(2)];
    sol0.v_k=[v0_k(3)];
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=[sol.u_k;u0_k(3)];
    v_k=[v0_k(1);v0_k(2);sol.v_k];

elseif n_helix+n_straight==3
    u_k = optimvar('u_k',2,1);
    v_k = optimvar('v_k',1,1);
    prob.Objective=helperfun([u_k(1);u_k(2);u0_k(3)],[v0_k(1);v0_k(2);v_k(1)]);
    sol0.u_k=[u0_k(1);u0_k(2)];
    sol0.v_k=[v0_k(3)];
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=[sol.u_k;u0_k(3)];
    v_k=[v0_k(1);v0_k(2);sol.v_k];

%     u_k = optimvar('u_k',3,1);
%     % v_k = optimvar('v_k',1,1);
%     prob.Objective=helperfun([u_k(1);u_k(2);u_k(3)],[v0_k(1);v0_k(2);v0_k(3)+0*u_k(3)]);
%     sol0.u_k=u0_k;
%     % sol0.v_k=[v0_k(3)];
%     sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
%     u_k=sol.u_k;
%     v_k=v0_k;
elseif n_helix+n_straight==2
    u_k = optimvar('u_k',2,1);
    prob.Objective=helperfun([u_k(1);u_k(2);u0_k(3)],[v0_k(1);v0_k(2);v0_k(3)]);
    sol0.u_k=[u0_k(1);u0_k(2)];
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=[sol.u_k;u0_k(3)];
    v_k=[v0_k(1);v0_k(2);v0_k(3)];
elseif n_helix+n_straight==1
    v_k = optimvar('v_k',1,1);
    prob.Objective=helperfun([u0_k(1);u0_k(2);u0_k(3)],[v0_k(1);v0_k(2);v_k(1)]);
    sol0.v_k=[v0_k(3)];
    sol = solve(prob,sol0,'Options',options,'Solver','fminunc');
    u_k=u0_k(3);
    v_k=[v0_k(1);v0_k(2);sol.v_k];
end


