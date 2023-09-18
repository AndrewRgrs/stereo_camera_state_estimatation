clc
clear all
%close all

load('dataset3.mat')


plotting = true; 

K1 = 1215;
K2 = 1715;
%K2 = 1300;

T_initial = zeros(4,4,K2-K1);
P_initial = zeros(6,6,K2-K1);


%Dead reckon to generate initial guess
for i = K1:K2 
    
    if i == K1
        %build transformation matrix using 6x1 ground truth column
        C_vk_i = vec2rot(theta_vk_i(:,K1))';
        T_initial(:,:,1) = [C_vk_i -C_vk_i*r_i_vk_i(:,K1); [0 0 0] 1];
        %Build initial vector and Rotation Matrix
        r(:,1) = r_i_vk_i(:,K1);
        C(:,:,1) = vec2rot(theta_vk_i(:,K1))';
        
        
    else
        time = t(i) - t(i-1);
        psi_k = w_vk_vk_i(:,i-1) * time;
        Psi = vec2rot(psi_k)';

        d_vk = v_vk_vk_i(:,i-1)* time;
        
        Xi = vec2tran([-d_vk;-psi_k]);
        
        C(:,:, i-K1 +1) = Psi* C(:,:, i-K1);
        r(:, i-K1 +1) = r(:, i-K1) + C(:,:, i-K1)'*d_vk;
        
        %Construct transformation matrix
        %T_initial(:,:,i-K1+1) = [ C(:,:, i-K1 +1), -C(:,:, i-K1 +1)*r(:, i-K1 +1); [0 0 0] 1];
        T_initial(:,:,i-K1+1) = Xi*T_initial(:,:,i-K1);
        C_vk_i = vec2rot(theta_vk_i(:,i));
        T_true(:,:,i-K1+1) = [C_vk_i -C_vk_i*r_i_vk_i(:,i); [0 0 0] 1];
        
        
        %Figure out uncertainty
        F = tranAd(T_initial(:,:,i-K1+1) \ (T_initial(:,:,i-K1)));
        Q_k= time*time * diag([v_var;w_var]);
        P_initial(:,:,i-K1+1) = F*P_initial(:,:,i-K1)*F' + Q_k; 
        
    end
    %4x4 transformation to 6x1 vector
    vecPose(:,i-K1+1) = tran2vec(T_initial(:,:,i-K1+1));
    %4x4 Transform to r and C
    x3(:,i-K1+1) = -(T_initial(1:3,1:3,i-K1+1)')*T_initial(1:3,4,i-K1+1);
    
    angle = rot2vec(C(:,:,i-K1 +1));
    vecPose2(:,i-K1 +1) = [r(:,i-K1 +1); angle];
end

%plot dead reckon stuff
if plotting
    x = vecPose(1:3,:);
    x2 = vecPose2(1:3,:);
    figure();
    plot3(rho_i_pj_i(1,:) , rho_i_pj_i(2,:), rho_i_pj_i(3,:), '.red')
    hold on
    plot3(x2(1,:),x2(2,:),x2(3,:) , 'm')
    hold on
    plot3(x3(1,:),x3(2,:),x3(3,:) ,'blue')
    plot3(r_i_vk_i(1,K1:K2),r_i_vk_i(2,K1:K2),r_i_vk_i(3,K1:K2), 'green');

    grid on
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    legend('','Transformation Model', 'Assignment Model', 'Ground Truth');
end
%% Plot How many Landmarks being seen

n_land = zeros(1 , (K2-K1)+1);
indg=0;
indr=0;

figure
hold on
for k = K1:K2
    for j = 1:20
        if y_k_j(1, k, j) ~= -1
            n_land(1,k-K1+1) = n_land(1,k-K1+1) +1;
        end
    end
    if n_land(k-K1+1) > 2
        plot(t(k),n_land(k-K1+1), '.green', 'markersize', 10);
        hold on

    else
        plot(t(k),n_land(k-K1+1), '.red', 'markersize', 10);
        hold on
    end
end
grid on
%plot(t_green(K1:K2),n_land);
xlabel('Time (s)')
ylabel('Number of visible landmarks M_k')
%% Setup Variables

T_op = T_initial;

P_o_prior = 0.1 * eye(6);

e_v_op = zeros(6, 1);
F = zeros(6, 6, K2 - K1);
Q = speye(6*(K2-K1),6*(K2-K1));
F_big = speye(6*(K2-K1),6*(K2-K1));

%e_y_op = zeros(4, 20, K2-K1);

%Z_j = zeros(4, 6, 20); %Combined measurement jacobian
%Z_k = zeros(80,6, K2-K1); %4*20 by 6 (jacobians for each landmark stacked)
Z_big = spalloc(80*(K2-K1), 6*(K2-K1), 80*6*(K2-K1)); 
S = zeros(4, 6); %Point to Camera frame mapping jacobian
G = zeros(4, 4); %Stereo Camera Jacobian
Z_k = [];

R = [];

%D = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
T_cv = [C_c_v -C_c_v*rho_v_c_v; [0 0 0 1]];
M = [fu 0 cu 0; 0 fv cv 0; fu 0 cu -fu*b; 0 fv cv 0];

R_j = diag (y_var);


input = -[v_vk_vk_i; w_vk_vk_i]; 
pert = ones(6*(K2-K1),1);

% form W matrix seperate so the while loop doesn't take forever
Qind =1;
for k = 1:(K2-K1)
    
    time = t(k+K1) - t(k+K1 -1);
    if k ==1
        Q(Qind:Qind+5,Qind:Qind+5) = P_o_prior;
    else
        Q(Qind:Qind+5,Qind:Qind+5) = time*time *diag([v_var; w_var]);
        Qind= Qind +6;
    end
        %rkind =1;
        R_k = [];
        for j= 1:20
            if y_k_j(1, k + K1, j) ~= -1 %landmark observable
                %R_k(rkind:rkind+3,rkind:rkind+3) = R_j;
                %rkind = rkind+4;
                R_k = blkdiag(R_k,R_j);
            else
                %rkind = rkind+4;
            end
        end
        R = sparse(blkdiag(R,R_k));
        %{
        if ~isempty(R_k)
            R(rind:rind+(length(R_k)-1), rind:rind+(length(R_k)-1)) = R_k;
            rind = rind + length(R_k);
            
        end
        %}
    
end
W = sparse(blkdiag(Q, R));
Winv = inv(W); %Assume this is gud
%% Do Gauss Newton
count = 0;
while (norm(pert) > 0.0001)
    Z_big = [];
    e_y_stacked = [];
    e_v_stacked = [];
    for k = 1: (K2-K1)
        if k ==1
            %e_v_op(:,1) = tran2vec(T_initial(:,:,1) / T_op(:,:,1));
            e_v_op = tran2vec(T_initial(:,:,1) / T_op(:,:,1));
            e_v_stacked = [e_v_stacked; e_v_op];
        else
            %Motion Model error section
            time = t(k+K1) - t(k+K1 -1);
            Xi = vec2tran(time*input(:,k+K1));
            %e_v_op(:,k) = tran2vec(Xi * T_op(:,:,k-1) / T_op(:,:,k) );
            e_v_op = tran2vec(Xi * T_op(:,:,k-1) / T_op(:,:,k) );
            e_v_stacked = [e_v_stacked; e_v_op];
            F(:,:,k-1) = tranAd( T_op(:,:,k) / T_op(:,:,k-1));

        end
        %Measurement Model Error Section
        %TODO: Remove Z_j for unobserable stuff
        Z_k = [];
        for j = 1:20 %iterate through all landmarks
            if y_k_j(1, k + K1, j) == -1 %landmark not observable
                
            else %landmark observable
                S = T_cv *point2fs( T_op(:,:,k) *[rho_i_pj_i(:,j);1] );
                G = camJac(M,T_cv*T_op(:,:,k)*[rho_i_pj_i(:,j);1]);
                %Z_j(:,:,j) = G *S;
                Z_k = [Z_k; G*S];
                
                %e_y_op(:, j, k) = y_k_j(:, k+K1, j) - cam(M, T_cv*T_op(:,:,k) *[rho_i_pj_i(:,j);1]);
                e_y_op = y_k_j(:, k+K1, j) - cam(M, T_cv*T_op(:,:,k) *[rho_i_pj_i(:,j);1]);
                e_y_stacked = [e_y_stacked;e_y_op];
                             
            end
            %Z_k(:,:,k) =reshape( permute(Z_j, [1 3 2]), [], 6);
            
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
    %z_ind1 = 1;
    %z_ind2 = 1;
    for i = 1:(K2-K1)

        F_big(f_ind:f_ind+5,f_ind:f_ind+5) = eye(6); 
        %Z_big(z_ind1:z_ind1 + 79,z_ind2:z_ind2+5) = Z_k(:,:,i); % Change this to so Z_big is 2x2
        %z_ind1 = z_ind1 +80;
        %z_ind2 = z_ind2 +6;
        if i ~= K2-K1
            F_big(f_ind+6:f_ind+11,f_ind:f_ind+5) = -F(:,:,i);
            f_ind = f_ind +6;

        end
    end
    H = [F_big; Z_big]; %should be 2x2
    %e_y_stacked = reshape(e_y_op, [], 1);
    %e_v_stacked = reshape(e_v_op, [], 1);
    e_big = [e_v_stacked; e_y_stacked];
    
    %L = chol(H' *Winv *H);
    %d = L\  H' *Winv *e_big;
    %pert = (L')\d;
    
    pert = (H'*Winv*H)\(H' *Winv *e_big);
    
    %update operating point
    ind =1;
    for k = 1: (K2-K1)
        T_change = vec2tran(pert(ind:ind+5));
        T_op(:,:,k) = T_change *T_op(:,:,k);
        ind = ind+6;
    end
    
    %force stop after 100 iterations
    count = count+1;
    %count = 100;
    if count >= 100
        pert = 0*pert;
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
covariance = inv(H' * Winv * H);
cov_ind = 1;
%extract (x y z) and C_estimated, calculate error and extract covariance
for i = K1:K2-1
    ind = i-K1+1; %index from 1 to (K2-K1)
    xyz(:,ind) = -(T_op(1:3,1:3,ind)')*T_op(1:3,4,ind);
    C_est(:,:,ind) = T_op(1:3,1:3,ind);
    
    C_true(:,:,ind) = vec2rot(theta_vk_i(:,i));
    
    %calculate error
    error_xyz(:,ind) = xyz(:,ind) - r_i_vk_i(:,i);
    error_C = eye(3) - C_est(:,:,ind) * C_true(:,:,ind);
    error_rot(:, ind) = [error_C(3,2); error_C(1,3); error_C(2,1)];
    
    %extract covariance for each timestep
    xyz_cov(:, ind) = diag(covariance(cov_ind:cov_ind+2,cov_ind:cov_ind+2));
    rot_cov(:, ind) = diag(covariance(cov_ind+3:cov_ind+5, cov_ind+3:cov_ind+5));
    cov_ind = cov_ind +6; %index for iterating through big 'ol covariance
   
end
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
axis([110 155 -0.17 0.17]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_xyz(2, :));
hold on
plot(t(K1:K2-2), 3*sqrt(xyz_cov(2,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(xyz_cov(2,1:end-1)), 'red');
legend
xlabel('Time [s]');
ylabel('Translational error in y [m]');
legend('Error $${\delta}r_{y,k}$$', '$$\pm 3{\sigma}_{r_{y,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.17 0.17]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_xyz(3, :));
hold on
plot(t(K1:K2-2), 3*sqrt(xyz_cov(3,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(xyz_cov(3,1:end-1)), 'red');
legend
xlabel('Time [s]');
ylabel('Translational error in z [m]');
legend('Error $${\delta}r_{z,k}$$', '$$\pm 3{\sigma}_{r_{z,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.17 0.17]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(1, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(1,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(1,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_x [rad]');
legend('Error $${\delta}{\theta}_{x,k}$$', '$$\pm 3{\sigma}_{{\theta}_{x,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.3 0.3]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(2, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(2,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(2,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_y [rad]');
legend('Error $${\delta}{\theta}_{y,k}$$', '$$\pm 3{\sigma}_{{\theta}_{y,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.3 0.3]);
figure('Units', 'centimeters','Position',[3 1 20 10])
plot(t(K1:K2-1), error_rot(3, :));
hold on
plot(t(K1:K2-2), 3*sqrt(rot_cov(3,1:end-1)), 'red');
plot(t(K1:K2-2), -3*sqrt(rot_cov(3,1:end-1)), 'red');
xlabel('Time [s]');
ylabel('Rotational error in \theta_z [rad]');
legend('Error $${\delta}{\theta}_{z,k}$$', '$$\pm 3{\sigma}_{{\theta}_{z,k}}$$', 'Interpreter', 'Latex')
axis([110 155 -0.3 0.3]);


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

