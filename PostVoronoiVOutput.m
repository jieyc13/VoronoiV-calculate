% Purpose: To calculate the Voronoi volume of each Voronoi cell of particles
%          Output the Voronoi volumns of particles in a *.dat file.
% Date: December 5th, 2018
% By Yucheng Jie
clear;
clc;

% Note: Need function MirrBnd.m and the width of mirror region 'delta' should be set properly in funcion MirrBnd.m

% Basic variables. 
% Grid: Nx*Ny*Nz. Domain size: Lx*Ly*Lz. Number of processes: PROC.
% Grid in each process: Nx*Ny*nz2. Number of particle: Np. Number of particle in each process: Np2.
% Type of particles: NtTotal. Select the particle 1:Nt to calculate the voronoi volumes. (Just set Nt=NtSphere for all types of particle.)
% Number of variables in 'particlefield*' file: NVar.
% Directory of data files: dir. Time steps for statistic: TimeSteps array.
Nx=384;
Ny=384;
Nz=384;
Lx=3.0;
Ly=1.5;
Lz=1;
PROC=191;
nz2=Nz/(PROC+1);
Np=6e5;
Np2=Np/(PROC+1);
Nt=8;
NtTotal=23;
NVar=13;
dir='/work/lihao/test/1200/part/data/';
TimeSteps=193000:1000:194000;

% Loop for every time step in array TimeSteps
for iFile=TimeSteps
    clearvars -except iFile Nx Ny Nz Lx Ly Lz PROC nz2 Np Np2 Nt NtTotal NVar dir;
	fprintf('iFile=%d\n',iFile);
	tic;
    % Variable columns: Variable [1,2,3] correspond to xp, yp and zp in 'particlefield*' files.
    VarClm=[1 2 3];
    % Load data in 'particlefield*' and save xp, yp, zp in matrix SumM.
    SumM=[];
    for i=0:PROC
        fdir=[dir,'particlefield',sprintf('%07d',iFile),'.',sprintf('%03d',i)];
        fid=fopen(fdir,'r');
        if (fid==-1)
            fprintf('File is not found.\n');
        end
        fseek(fid,4,-1);                                    % Useing matlab to read binary data files, there is a difference in the front of files.
        [Matrix]=fread(fid,'real*8');
        Matrix=reshape(Matrix,Np2,NtTotal,NVar);
        SumM(1+i*Np2:(i+1)*Np2,:,:)=Matrix(:,:,VarClm);     % Save xp, yp, zp in matrix SumM.
        fclose(fid);
    end
    
    Par=cell(NtTotal,1);
    xp=cell(NtTotal,1);
    yp=cell(NtTotal,1);
    zp=cell(NtTotal,1);
    
    % Assignment for xp, yp, zp of different types of particles.
    for i=1:NtTotal
        Par{i}=squeeze(SumM(:,i,:));
    end
    for i=1:NtTotal
        xp{i}=Par{i}(:,1);
        yp{i}=Par{i}(:,2);
        zp{i}=Par{i}(:,3);
    end

    % Select particle type 1:Nt to calculate the voronoi volume for every type.
    volume=cell(Nt,1);
    VNorm=cell(Nt,1);
    for i=1:Nt
        fprintf('Particle type:%d\n',i);
        x=xp{i};
        y=yp{i};
        z=zp{i};

        xyz=[x, y, z];

        % Mirrored boundary for particles. Note that the width of mirror region 'delta' should be set properly in funcion MirrBnd.m
        Pos=MirrBnd(xyz,Lx,Ly,Lz);

        % Voronoi diagram for each particle
        [V,C] = voronoin(Pos);%,{'Qbb'});

        % Index of original particles.
        iIn=1:1:Np;

        % Calculate the voronoi volume of each particle
        volume{i}=zeros(size(iIn));
        for ip=iIn(1):iIn(end)
            % Use convhull function to calculate the voronoi volumns.
            [K,volume{i}(ip)]=convhull(V(C{ip},1),V(C{ip},2),V(C{ip},3));
        end 
        % Normalized the Voronoi volume so that V=1 stands for the size of Voronoi cell in inertial random distribution of particles.
        VNorm{i}=volume{i}*Np/(Lx*Ly*Lz);
    end

    %% Output Voronoi volume of each particle
    for iNt=1:Nt
        fid=fopen(['VoronoiV_',sprintf('%07d',iFile),'_',sprintf('%02d',iNt),'.dat'],'w');
        for i=iIn(1):iIn(end)
            fprintf(fid,'%e %e %e %e\n',xp{iNt}(i),yp{iNt}(i),zp{iNt}(i),VNorm{iNt}(i));
        end
        fclose(fid);
    end
	toc;
end

    













