function [ A ] = MirrBnd( pos, Lx, Ly, Lz )
%MirrBnd Mirrored boundary of the particle position
%   Pre-post for voronoi 

% Note: delta is the width of mirror boundary region, while should be set properly at first.
       % Origal domain is [0,Lx]*[0,Ly]*[0,*Lz]. (Must start from 0. The function is designed for this domain.)
       % Input: pos-[x,y,z] matrix of particle location.
       % Output: A-[-delta, Lx+delta]*[-delta, Ly+delta]*[-delta, Lz+delta]. 
       %        Using mirror boundary condition to extend the original domain.
delta=1e-1;
A=pos;

% Cycle in all direction to extend original domain.
for ix=-1:1
    for iy=-1:1
        for iz=-1:1
            if ((ix==0)&&(iy==0)&&(iz==0)) 
                continue;
            end
            % Find proper paritlces near the specific boundary.
            xAdd=find( ((-ix*pos(:,1)>=-(ix+1.)/2*Lx) & (-ix*pos(:,1)<=-(ix+1.)/2*Lx + delta)) | ix==0 );
            yAdd=find( ((-iy*pos(:,2)>=-(iy+1.)/2*Ly) & (-iy*pos(:,2)<=-(iy+1.)/2*Ly + delta)) | iy==0 );
            zAdd=find( ((-iz*pos(:,3)>=-(iz+1.)/2*Lz) & (-iz*pos(:,3)<=-(iz+1.)/2*Lz + delta)) | iz==0 );
            iAdd=intersect(intersect(xAdd,yAdd),zAdd);      % Intersection
            % Using mirror boundary condition to extend the original domain.
            A=[A; ...
                [(abs(2*ix+1)-1)*Lx - 2*(abs(ix)-1./2)*pos(iAdd,1),...
                 (abs(2*iy+1)-1)*Ly - 2*(abs(iy)-1./2)*pos(iAdd,2),...
                 (abs(2*iz+1)-1)*Lz - 2*(abs(iz)-1./2)*pos(iAdd,3)] ];
        end
    end
end

end

