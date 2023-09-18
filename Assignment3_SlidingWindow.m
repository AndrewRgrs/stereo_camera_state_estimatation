clc
clear all
close all

load('dataset3.mat')
% Time to do shit
M = [fu 0 cu 0; 0 fv cv 0; fu 0 cu -fu*b; 0 fv cv 0];

plotting = true; 

K1 = 1215;
K2 = 1715;
kappa = 50;

T_est = zeros(4,4,(K2-K1));
Cov_est = zeros(6,6,(K2-K1));

for k = K1:K2-1
    if k == K1
        C_vk_i = vec2rot(theta_vk_i(:,k))';
        T_start = [C_vk_i -C_vk_i*r_i_vk_i(:,k); [0 0 0] 1];

        [T_initial] = deadreckon(T_start,k,k + kappa);
        [T_est(:,:,1),Cov_est(:,:,1)] = GN(T_initial,zeros(6), k, k+kappa, M);
    else
        T_start = T_est(:,:, k-K1); %Previous Estimated State
        P_start = Cov_est(:,:, k-K1); %Covariance of previous estimate
        
        T_initial = deadreckon(T_start,k, k+kappa);
        [T_est(:,:,k-K1+1),Cov_est(:,:,k-K1+1)] = GN(T_initial, P_start, k, k+kappa, M);
        
    end
end


%% plot trajectory and error
xyz = zeros(3,(K2-K1));
C_est = zeros(3,3,(K2-K1));
C_true = zeros(3,3, (K2-K1));
error_xyz = zeros(3,(K2-K1));
error_rot = zeros(3,(K2-K1));
xyz_cov = zeros(3,(K2-K1));
rot_cov = zeros(3,(K2-K1));

%Find Covariance Matirx
%extract (x y z) and C_estimated, calculate error and extract covariance
for i = K1:K2-1
    ind = i-K1+1; %index from 1 to (K2-K1)
    xyz(:,ind) = -(T_est(1:3,1:3,ind)')*T_est(1:3,4,ind);
    C_est(:,:,ind) = T_est(1:3,1:3,ind);
    
    C_true(:,:,ind) = vec2rot(theta_vk_i(:,i));
    
    %calculate error
    error_xyz(:,ind) = xyz(:,ind) - r_i_vk_i(:,i);
    error_C = eye(3) - C_est(:,:,ind) * C_true(:,:,ind);
    error_rot(:, ind) = [error_C(3,2); error_C(1,3); error_C(2,1)];
    
    %extract covariance for each timestep
    xyz_cov(:, ind) = diag(Cov_est(1:3,1:3,ind));
    rot_cov(:, ind) = diag(Cov_est(4:6,4:6,ind));
   
end
%% Plot
close all
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_xyz(1, :));
hold on
plot(t(K1:K2-2), 3*sqrt(xyz_cov(1,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(xyz_cov(1,1:end-1)), 'red');
legend
xlabel('Time [s]');
ylabel('Translational error in x [m]');
legend('Error $${\delta}r_{x,k}$$', '$$\pm 3{\sigma}_{r_{x,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -1 1]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_xyz(2, :));
hold on
plot(t(K1:K2-2), 3*sqrt(xyz_cov(2,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(xyz_cov(2,1:end-1)), 'red');
legend
xlabel('Time [s]');
ylabel('Translational error in y [m]');
legend('Error $${\delta}r_{y,k}$$', '$$\pm 3{\sigma}_{r_{y,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -1 1]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_xyz(3, :));
hold on
plot(t(K1:K2-2), 3*sqrt(xyz_cov(3,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(xyz_cov(3,1:end-1)), 'red');
legend
xlabel('Time [s]');
ylabel('Translational error in z [m]');
legend('Error $${\delta}r_{z,k}$$', '$$\pm 3{\sigma}_{r_{z,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -1 1]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(1, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(1,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(1,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_x [rad]');
legend('Error $${\delta}{\theta}_{x,k}$$', '$$\pm 3{\sigma}_{{\theta}_{x,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.-1 1]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(2, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(2,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(2,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_y [rad]');
legend('Error $${\delta}{\theta}_{y,k}$$', '$$\pm 3{\sigma}_{{\theta}_{y,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.-1 1]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(3, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(3,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(3,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_z [rad]');
legend('Error $${\delta}{\theta}_{z,k}$$', '$$\pm 3{\sigma}_{{\theta}_{z,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.-1 1]);


figure
plot3(rho_i_pj_i(1,:) , rho_i_pj_i(2,:), rho_i_pj_i(3,:), '.red')
hold on
plot3(xyz(1,:),xyz(2,:),xyz(3,:) ,'blue')
plot3(r_i_vk_i(1,K1:K2),r_i_vk_i(2,K1:K2),r_i_vk_i(3,K1:K2), 'green');
grid on

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Landmarks','GN Estimate', 'Ground Truth','Location','northeast')

%% functions


function [T_initial] = deadreckon(T_start, k1, k2)
    load dataset3.mat ;
    T_initial = zeros(4,4,k2-k1);
   

    %Dead reckon to generate initial guess
    for i = k1:k2
        if i == k1
            T_initial(:,:,1) = T_start;
            %Build initial vector and Rotation Matrix


        else
            time = t(i) - t(i-1);
            psi_k = w_vk_vk_i(:,i-1) * time;

            d_vk = v_vk_vk_i(:,i-1)* time;

            Xi = vec2tran([-d_vk;-psi_k]);

            %Construct transformation matrix
            %T_initial(:,:,i-K1+1) = [ C(:,:, i-K1 +1), -C(:,:, i-K1 +1)*r(:, i-K1 +1); [0 0 0] 1];
            T_initial(:,:,i-k1+1) = Xi*T_initial(:,:,i-k1);

        end

    end

end

function [T_k1, Cov] = GN(T_initial, P, k1, k2, M)
    load dataset3.mat;
    T_op = T_initial;
    P_o_prior = P;

    F = zeros(6, 6, k2 - k1);
    Q = speye(6*(k2-k1),6*(k2-k1));
    F_big = speye(6*(k2-k1),6*(k2-k1));

    Z_big = spalloc(80*(k2-k1), 6*(k2-k1), 80*6*(k2-k1)); 
    S = zeros(4, 6); %Point to Camera frame mapping jacobian
    G = zeros(4, 4); %Stereo Camera Jacobian
    Z_k = [];

    R = [];

    T_cv = [C_c_v -C_c_v*rho_v_c_v; [0 0 0 1]];
    

    R_j = diag (y_var);

    input = -[v_vk_vk_i; w_vk_vk_i]; 
    pert = ones(6*(k2-k1),1);

    Qind =1;
%Make W matrix
    for k = 1:(k2-k1)
        time = t(k+k1) - t(k+k1 -1);
        if k ==1
            Q(Qind:Qind+5,Qind:Qind+5) = P_o_prior;
        else
            Q(Qind:Qind+5,Qind:Qind+5) = time*time *diag([v_var; w_var]);
            Qind= Qind +6;
        end
            %rkind =1;
            R_k = [];
            for j= 1:20
                if y_k_j(1, k + k1, j) ~= -1 %landmark observable

                    R_k = blkdiag(R_k,R_j);
                end
            end
            R = sparse(blkdiag(R,R_k));

    end
    W = sparse(blkdiag(Q, R));
    Winv = inv(W); %Assume this is gud
    % Do Gauss Newton
    count = 0;
    while (norm(pert) > 0.001)
        Z_big = [];
        e_y_stacked = [];
        e_v_stacked = [];
        for k = 1: (k2-k1)
            if k ==1
                e_v_op = tran2vec(T_initial(:,:,1) / T_op(:,:,1));
                e_v_stacked = [e_v_stacked; e_v_op];
            else
                %Motion Model error section
                time = t(k+k1) - t(k+k1 -1);
                Xi = vec2tran(time*input(:,k+k1));

                e_v_op = tran2vec(Xi * T_op(:,:,k-1) / T_op(:,:,k) );
                e_v_stacked = [e_v_stacked; e_v_op];
                F(:,:,k-1) = tranAd( T_op(:,:,k) / T_op(:,:,k-1));

            end
            %Measurement Model Error Section
            Z_k = [];
            for j = 1:20 %iterate through all landmarks
                if y_k_j(1, k + k1, j) == -1 %landmark not observable

                else %landmark observable
                    S = T_cv *point2fs( T_op(:,:,k) *[rho_i_pj_i(:,j);1] );
                    G = camJac(M,T_cv*T_op(:,:,k)*[rho_i_pj_i(:,j);1]);
                    Z_k = [Z_k; G*S];

                    e_y_op = y_k_j(:, k+k1, j) - cam(M, T_cv*T_op(:,:,k) *[rho_i_pj_i(:,j);1]);
                    e_y_stacked = [e_y_stacked;e_y_op];

                end

            end
            if ~isempty(Z_k)
                Z_big = sparse(blkdiag(Z_big,Z_k));
            else
                Z_big = sparse([Z_big, zeros( size(Z_big,1),6)]); %append columns of zeros 
            end
        end
        %Build top of H matirx - F_big
        %Inefficient method that works:
        f_ind = 1;

        for i = 1:(k2-k1)

            F_big(f_ind:f_ind+5,f_ind:f_ind+5) = eye(6); 

            if i ~= k2-k1
                F_big(f_ind+6:f_ind+11,f_ind:f_ind+5) = -F(:,:,i);
                f_ind = f_ind +6;

            end
        end
        H = [F_big; Z_big]; %should be 2x2
        e_big = [e_v_stacked; e_y_stacked];

        pert = (H'*Winv*H)\(H' *Winv *e_big);

        %update operating point
        ind =1;
        for k = 1: (k2-k1)
            T_change = vec2tran(pert(ind:ind+5));
            T_op(:,:,k) = real(T_change *T_op(:,:,k));
            ind = ind+6;
        end

        %force stop after 100 iterations
        count = count+1;
        %count = 100;
        if count >= 100
            pert = 0*pert;
        end

    end
    T_k1 = real(T_op(:,:,1));
    covariance = inv(H'*Winv*H);
    Cov = covariance(1:6,1:6);

end