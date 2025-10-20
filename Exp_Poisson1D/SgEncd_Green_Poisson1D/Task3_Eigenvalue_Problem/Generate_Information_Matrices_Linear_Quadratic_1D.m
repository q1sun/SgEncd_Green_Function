function [P_mesh,T_mesh,V_basis,T_basis] = Generate_Information_Matrices_Linear_Quadratic_1D(left,right,h,basis_type)

format short e

%% 1. generate information matrices P and T for a partition of [left,right]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------%
% 1-1. linear elements (also the grid points)
%-----------------------------------------------------------%
N_linear = (right-left)/h;
P_linear = zeros(1,N_linear+1);
T_linear = zeros(2,N_linear);

for i = 1 : N_linear+1
    P_linear(1,i) = left + (i-1)*h;
end

for i = 1 : N_linear
    T_linear(1,i) = i;
    T_linear(2,i) = i+1;
end

% boundary nodes (global index among all finite element nodes)
BN_linear = [1,length(P_linear)];
%-----------------------------------------------------------%
% 1-2. quadratic elements
%-----------------------------------------------------------%
N_quadratic = (right-left)/h;
dh = h/2;
dN = N_quadratic*2;
P_quadratic = zeros(1,dN+1);
T_quadratic = zeros(2,dN);

for i = 1 : dN+1
    P_quadratic(1,i) = left+(i-1)*dh;
end

for i = 1 : N_quadratic
    T_quadratic(1,i) = 2*i-1;
    T_quadratic(2,i) = 2*i+1;
    T_quadratic(3,i) = 2*i;
end

% boundary nodes (global index among all finite element nodes)
BN_quadratic = [1,length(P_quadratic)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 2. display mesh and information matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if plot_mesh == 1
%     
%     figure
%     
%     % linear elements
%     subplot(2,1,1)
%     V_linear = [P_linear; zeros(1,size(P_linear,2))];    
%     patch('Faces',T_linear','Vertices',V_linear','FaceColor','w');
%     hold on
%     plot(V_linear(1,:),V_linear(2,:),'*');
%     text(V_linear(1,:),V_linear(2,:)-0.05,num2str((1:size(P_linear,2))'),'FontSize',12,'Color','red');
%     % order of face displayed in the mean value coordinates of the element    
%     V_text = [mean([P_linear(T_linear(1,:));P_linear(T_linear(2,:))]); ones(1,size(T_linear,2))*0.1 ];    
%     text(V_text(1,:),V_text(2,:),num2str((1:size(T_linear,2))'),'FontSize',10,'HorizontalAlignment','left','BackgroundColor',[.7 .9 .7]);
%     % some parameters
%     axis([min(P_linear)-0.02 max(P_linear)+0.02 -0.5 0.5]);
%     title('information matrices for linear elements (also the mesh matrices)');
%     xlabel('x label');
%     ylabel('y label');
%     hold off
%     
%     % quadratic elements
%     subplot(2,1,2)
%     V_quadratic = [P_quadratic; zeros(1,size(P_quadratic,2))];    
%     patch('Faces',T_linear','Vertices',V_linear','FaceColor','w');
%     hold on    
%     plot(V_quadratic(1,:),V_quadratic(2,:),'*');
%     text(V_quadratic(1,:),V_quadratic(2,:)-0.05,num2str((1:size(P_quadratic,2))'),'FontSize',12,'Color','red');
%     % order of face displayed in the mean value coordinates of the element    
%     V_text = [mean([P_linear(T_linear(1,:));P_linear(T_linear(2,:))]); ones(1,size(T_linear,2))*0.1 ];    
%     text(V_text(1,:),V_text(2,:),num2str((1:size(T_linear,2))'),'FontSize',10,'HorizontalAlignment','left','BackgroundColor',[.7 .9 .7]);
%     % some parameters
%     axis([min(P_linear)-0.02 max(P_linear)+0.02 -0.5 0.5]);
%     title('information matrices for quadratic elements');
%     xlabel('x label');
%     ylabel('y label');
%     hold off
%             
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh matrices are the same as the linear elements
P_mesh = P_linear;
T_mesh = T_linear;
% finite elements
if basis_type == 101 % linear elements
    V_basis = P_linear;
    T_basis = T_linear;
elseif basis_type == 102 % quadratic elements
    V_basis = P_quadratic;
    T_basis = T_quadratic; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
