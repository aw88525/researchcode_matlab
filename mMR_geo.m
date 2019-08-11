%This code is to generate the 3D coordinate for the voxels and crystal
%positions for ray tracing algorithm
%Scanner type is cylinder scanner
%Rsector displacement dx,dy,dz
%Rsector # n_rsec
clear;clc;close all;
disp_x = 0;
disp_y = 300+0; %mm
disp_z = 0;
%first gap position
n_rsec = 56;
N_mx = 1;
N_mz = 8;  %used to be 4
lm_x = 0;
lm_z = 32; %mm
N_cx = 8;
N_cz = 8;
lc_x = 4; %mm
lc_y = 20; %mm
lc_z = 4; %mm

%an array for storing the the center points of the rsectors
c_rsec = zeros(3,n_rsec);
%to compute the center points
for i=1:n_rsec
    c_rsec(1,i) = disp_x*cos(2*pi*(i-1)/n_rsec) + disp_y*sin(2*pi*(i-1)/n_rsec);
    c_rsec(2,i) = -disp_x*sin(2*pi*(i-1)/n_rsec) + disp_y*cos(2*pi*(i-1)/n_rsec);
    c_rsec(3,i) = 0;
end
% plot(c_rsec(1,:),c_rsec(2,:));


figure;
scatter3(c_rsec(1,:),c_rsec(2,:),c_rsec(3,:), 'g.');
% hold on;
% scatter3(c_gap(1,:),c_gap(2,:),c_gap(3,:), 'b.');
%Now calculate the center points for each module
%The repeat number for the modules in the x coordinate denoted as N_mx
%The repeat number for the modules in the y coordinate denoted as N_my
%The repeat number for the modules in the y coordinate denoted as N_mz
if disp_x ~= 0 
    N_mx = 0;
elseif disp_y ~= 0 
    N_my = 0;
elseif disp_z  ~= 0
    N_mz = 0;
else
    error('error occurred!');
end

%l_mx,l_my, l_mz for repeating intervals for each module
%in this case N_my = 0, so let us just consider x and z coordinate
c_module = zeros(3,N_mx*N_mz*n_rsec);
if(mod(N_mx,2) == 0 && mod(N_mz,2) == 0)    
   for k=1:n_rsec
    for i=-N_mx/2+1:N_mx/2
        for j=-N_mz/2+1:N_mz/2
             c_module(1,j+N_mz/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(1,k) - (0.5*lm_x + (i-1)*lm_x)*cos(2*pi*(k-1)/n_rsec);
             c_module(2,j+N_mz/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(2,k) + (0.5*lm_x + (i-1)*lm_x)*sin(2*pi*(k-1)/n_rsec);
             c_module(3,j+N_mz/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(3,k) - 0.5*lm_z - (j-1)*lm_z;
        end
    end
   end
elseif(mod(N_mx,2) == 1 && mod(N_mz,2) == 0)
   for k=1:n_rsec
    for i=-(N_mx-1)/2:(N_mx-1)/2
        for j=-N_mz/2+1:N_mz/2
             c_module(1,j+N_mz/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(1,k) - (i*lm_x)*cos(2*pi*(k-1)/n_rsec);
             c_module(2,j+N_mz/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(2,k) + (i*lm_x)*sin(2*pi*(k-1)/n_rsec);
             c_module(3,j+N_mz/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(3,k) - 0.5*lm_z - (j-1)*lm_z;
        end
    end
   end 
elseif(mod(N_mx,2) == 0 && mod(N_mz,2) == 1)
    for k=1:n_rsec
    for i=-N_mx/2+1:N_mx/2
        for j=-(N_mz-1)/2:(N_mz-1)/2
             c_module(1,j+(N_mz+1)/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(1,k) - (0.5*lm_x + (i-1)*lm_x)*cos(2*pi*(k-1)/n_rsec);
             c_module(2,j+(N_mz+1)/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(2,k) + (0.5*lm_x + (i-1)*lm_x)*sin(2*pi*(k-1)/n_rsec);
             c_module(3,j+(N_mz+1)/2+(i+N_mx/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(3,k) -j*lm_z;
        end
    end
    end 
else
   for k=1:n_rsec
    for i=-(N_mx-1)/2:(N_mx-1)/2
        for j=-(N_mz-1)/2:(N_mz-1)/2
             c_module(1,j+(N_mz+1)/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(1,k) - (i*lm_x)*cos(2*pi*(k-1)/n_rsec);
             c_module(2,j+(N_mz+1)/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(2,k) + (i*lm_x)*sin(2*pi*(k-1)/n_rsec);
             c_module(3,j+(N_mz+1)/2+(i+(N_mx+1)/2-1)*N_mz+(k-1)*N_mx*N_mz) = c_rsec(3,k) - j*lm_z;
        end
    end
   end   
end
    
    
    
%finally calculate the center points of each crystal
%l_cx,l_cy, l_cz for repeating intervals for each module
%in this case N_cy = 0, so let us just consider x and z coordinate
c_crystal = zeros(3,N_cx*N_cz*N_mx*N_mz*n_rsec);
if(mod(N_cx,2) == 0 && mod(N_cz,2) == 0)    
   for p=1:n_rsec
       for q=1:N_mx*N_mz 
        for i=-N_cx/2+1:N_cx/2
            for j=-N_cz/2+1:N_cz/2
                 c_crystal(1,j+N_cz/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(1,q+(p-1)*N_mx*N_mz) - (0.5*lc_x + (i-1)*(lc_x))*cos(2*pi*(p-1)/n_rsec);
                 c_crystal(2,j+N_cz/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(2,q+(p-1)*N_mx*N_mz) + (0.5*lc_x + (i-1)*(lc_x))*sin(2*pi*(p-1)/n_rsec);
                 c_crystal(3,j+N_cz/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(3,q+(p-1)*N_mx*N_mz) - 0.5*lc_z - (j-1)*lc_z;
            end
            
        end
%         c_crystal(1,j+N_cz/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz+gap_num) = 11.0922*cos(2*pi*(p-1)/n_rsec) + 84.2536*sin(2*pi*(p-1)/n_rsec);
%         gap_num = gap_num+1;
       end
   end
elseif(mod(N_cx,2) == 1 && mod(N_cz,2) == 0)
   for p=1:n_rsec
       for q=1:N_mx*N_mz 
        for i=-(N_cx-1)/2:(N_cx-1)/2
            for j=-N_cz/2+1:N_cz/2
                 c_crystal(1,j+N_cz/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(1,q+(p-1)*N_mx*N_mz) - (i*(lc_x))*cos(2*pi*(p-1)/n_rsec);
                 c_crystal(2,j+N_cz/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(2,q+(p-1)*N_mx*N_mz) + (i*(lc_x))*sin(2*pi*(p-1)/n_rsec);
                 c_crystal(3,j+N_cz/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(3,q+(p-1)*N_mx*N_mz) - 0.5*lc_z - (j-1)*lc_z;
            end
           
        end
       end 
   end
elseif(mod(N_cx,2) == 0 && mod(N_cz,2) == 1)
    for p=1:n_rsec
       for q=1:N_mx*N_mz 
        for i=-N_cx/2+1:N_cx/2
            for j=-(N_cz-1)/2:(N_cz-1)/2
                 c_crystal(1,j+(N_cz+1)/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(1,q+(p-1)*N_mx*N_mz) - (0.5*lc_x + (i-1)*(lc_x))*cos(2*pi*(p-1)/n_rsec);
                 c_crystal(2,j+(N_cz+1)/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(2,q+(p-1)*N_mx*N_mz) + (0.5*lc_x + (i-1)*(lc_x))*sin(2*pi*(p-1)/n_rsec);
                 c_crystal(3,j+(N_cz+1)/2+(i+N_cx/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(3,q+(p-1)*N_mx*N_mz) -j*lc_z;
            end
            
        end
       end 
    end
else
   for p=1:n_rsec
       for q=1:N_mx*N_mz 
        for i=-(N_cx-1)/2:(N_cx-1)/2
            for j=-(N_cz-1)/2:(N_cz-1)/2
                 c_crystal(1,j+(N_cz+1)/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(1,q+(p-1)*N_mx*N_mz) - (i*(lc_x))*cos(2*pi*(p-1)/n_rsec);
                 c_crystal(2,j+(N_cz+1)/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(2,q+(p-1)*N_mx*N_mz) + (i*(lc_x))*sin(2*pi*(p-1)/n_rsec);
                 c_crystal(3,j+(N_cz+1)/2+(i+(N_cx+1)/2-1)*N_cz+(q-1)*N_cx*N_cz+(p-1)*N_mx*N_mz*N_cx*N_cz) = c_module(3,q+(p-1)*N_mx*N_mz) - j*lc_z;
            end
           
        end
       end
   end
end
 


zmin = min(c_crystal(3,:));
x = ones(1,size(c_crystal,2));
c_crystal_zmin = c_crystal(:,x & c_crystal(3,:) == zmin); 


% c_crystal_center_zmin(:,1:8) = c_crystal_zmin(:,1:8);
% c_crystal_center_zmin(2,1:8) = c_crystal_center_zmin(2,1:8)+7;
% 
% 
% for i=2:n_rsec
%   for j=1:8
%     c_crystal_center_zmin(1,(i-1)*8+j) = c_crystal_center_zmin(1,j)*cos(2*pi*(i-1)/n_rsec) + c_crystal_center_zmin(2,j)*sin(2*pi*(i-1)/n_rsec);
%     c_crystal_center_zmin(2,(i-1)*8+j) = -c_crystal_center_zmin(1,j)*sin(2*pi*(i-1)/n_rsec) + c_crystal_center_zmin(2,j)*cos(2*pi*(i-1)/n_rsec);
%     c_crystal_center_zmin(3,(i-1)*8+j) = 0;
%   end
% end
% 
% 
% 
% 
% 
% for i=9:1:192
%    if(mod(i,8)~=0)
%      c_crystal_center_zmin(:,i) = c_crystal_center_zmin(:,mod(i,8))+c_rsec(:,floor((i-1)/8)+1)-c_rsec(:,1);
%    else
%      c_crystal_center_zmin(:,i) = c_crystal_center_zmin(:,8)+c_rsec(:,floor((i-1)/8)+1)-c_rsec(:,1);
%    end
% end



index_array = 1:448;
inversed_index_array=[];
for i=1:8:448
    temp = fliplr(i:i+7);
    inversed_index_array = [inversed_index_array temp];
end
shift_inversed_index_array(1:444) = inversed_index_array(5:448);
shift_inversed_index_array(445:448) = inversed_index_array(1:4); %This stores the indices for the crystal in one ring
% shift_inversed_index_array = inversed_index_array;

z_crystal = sort(c_crystal(3,:));



for i=1:length(z_crystal)
    c_crystal_new(3,i) = z_crystal(i);
    c_crystal_new(1,i) = c_crystal_zmin(1,shift_inversed_index_array(mod(i-1,448)+1));
    c_crystal_new(2,i) = c_crystal_zmin(2,shift_inversed_index_array(mod(i-1,448)+1));
%     c_crystal_newnew(3,i) = z_crystal(i);
%     c_crystal_newnew(1,i) = c_crystal_center_zmin(1,shift_inversed_index_array(mod(i-1,192)+1));
%     c_crystal_newnew(2,i) = c_crystal_center_zmin(2,shift_inversed_index_array(mod(i-1,192)+1));
%     c_crystal_new(4,i) = fix((shift_inversed_index_array(mod(i-1,192)+1)-1)/8);
end
% 



% c_crystal_new = zeros(3, size(z_crystal_new,2));


% 



% c_crystal_new_gap(3,1:192) = (c_crystal_new(3,577:768)+c_crystal_new(3,769:960))./2;
% c_crystal_new_gap(1,1:192) = c_crystal_new(1,577:768);
% c_crystal_new_gap(2,1:192) = c_crystal_new(2,577:768);
% 
% c_crystal_new_gap(3,193:384) = double(int32(c_crystal_new(3,1345:1536)+c_crystal_new(3,1537:1728))./2);
% c_crystal_new_gap(1,193:384) = c_crystal_new(1,577:768);
% c_crystal_new_gap(2,193:384) = c_crystal_new(2,577:768);
% 
% c_crystal_new_gap(3,385:576) =  (c_crystal_new(3,2113:2304)+c_crystal_new(3,2305:2496))./2;
% c_crystal_new_gap(1,385:576) = c_crystal_new(1,577:768);
% c_crystal_new_gap(2,385:576) = c_crystal_new(2,577:768);

fid_x = fopen('C:\Users\wsy88\Documents\MATLAB\crystal_position_mMR_surface_points.txt','wt');
fprintf(fid_x,'%f\t%f\t%f\n',c_crystal_new);
% fid_xx = fopen('C:\Users\tech\Documents\MATLAB\crystal_position_VersaPET_fourringsnogaps_center_shift.txt','wt');
% fprintf(fid_xx,'%f\t%f\t%f\n',c_crystal_newnew);

figure;
% scatter3(c_module(1,:),c_module(2,:),c_module(3,:), 'g.');
% hold on;
scatter3(c_crystal_new(1,:),c_crystal_new(2,:),c_crystal_new(3,:),'r.');

hold on;

fid_xx = fopen('C:\Users\wsy88\Documents\MATLAB\tof_center.txt','r');
tof_center = fscanf(fid_xx,'%f\t%f\t%f\n');
tof_center = reshape(tof_center', [3, 29295/3])';

scatter3(tof_center(:,1), tof_center(:,2), tof_center(:,3), 'g.');





%Now compute the coordinate for the voxel grid 
%It is best to create a voxel grid with the boundary as the detectors
% 
% limit_x = fix(disp_y/voxel_l)+1;
% limit_y = limit_x;
% limit_z = fix(((1.5*lm_z+lc_z*2)/voxel_l))+1;
% 
% voxel_x = (-limit_x:limit_x).*voxel_l;
% voxel_y = (-limit_y:limit_y).*voxel_l;
% voxel_z = (-limit_z:limit_z).*voxel_l;
% N_vx = size(voxel_x,2);
% N_vy = size(voxel_y,2);
% N_vz = size(voxel_z,2);
% voxel_pos = zeros(3,N_vx*N_vy*N_vz);
% 
% for i=1:N_vx
%     for j=1:N_vy
%         for k=1:N_vz
%             voxel_pos(1,k+(j-1)*N_vz+(i-1)*N_vy*N_vz) = voxel_x(i);
%             voxel_pos(2,k+(j-1)*N_vz+(i-1)*N_vy*N_vz) = voxel_y(j);
%             voxel_pos(3,k+(j-1)*N_vz+(i-1)*N_vy*N_vz) = voxel_z(k);
%         end
%     end
% end
% 
% 
% fid = fopen('C:\Users\tech\Documents\MATLAB\crystal_position_center.txt','wt');
% fprintf(fid,'%f\t%f\t%f\n',c_crystal_new);
% % c_crystal_new(3,1:192) = -20.8125;
% fid2 = fopen('C:\Users\tech\Documents\MATLAB\crystal_position_xyz_surface.txt','wt');
% fprintf(fid2,'%f\t%f\t\%f\n',c_crystal_new(:,1:192));



