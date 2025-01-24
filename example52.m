%% example 5.2
% comparisson for VEM and SSE-VEM applied to hexagonal and voronoi meshes
% exact solution and rhs
uex = @(x,y) sin(pi*x).*sin(pi*y);
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);

%% SSE-VEM Voronoi
% read meshes and number of files
files         = dir('./meshes/voro*');
numFiles      = size(files,1); 
errSSEVEMvoro = nan(numFiles,1);
hhSSEVEMvoro  = nan(numFiles,1);
errVEMvoro    = nan(numFiles,1);
hhVEMvoro     = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % exact solution at nodes
    ue = uex(mesh.verts(:,1),mesh.verts(:,2));
    % assemble matrix and rhs for VEM and SSE-VEM
    [AVEM,ASSEVEM,b,out,A2] = SSEVEM(mesh,f);
    % SSE-VEM solution
    freeNode = setdiff(1:size(mesh.verts,1),mesh.bdNodes)';
    u = zeros(size(mesh.verts,1),1);
    u(freeNode) = ASSEVEM(freeNode,freeNode)\b(freeNode);
    % compute L2-error
    errSSEVEMvoro(pt+1) = getL2Error(out,uex,u); % (u-ue)'*A2*(u-ue);%
    % VEM solution
    u = zeros(size(mesh.verts,1),1);
    u(freeNode) = AVEM(freeNode,freeNode)\b(freeNode);
    % compute L2-error
    errVEMvoro(pt+1) = getL2Error(out,uex,u); % (u-ue)'*A2*(u-ue);%
    %  compute diam
    hhSSEVEMvoro(pt+1) = mesh.hmin;
    hhVEMvoro(pt+1) = mesh.hmin;
    pt = pt + 1;
end

%% SSE-VEM Hexagons
% read meshes and number of files
files         = dir('./meshes/hexa*');
numFiles      = size(files,1); 
errSSEVEMhexa = nan(numFiles,1);
hhSSEVEMhexa  = nan(numFiles,1);
errVEMhexa    = nan(numFiles,1);
hhVEMhexa     = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % exact solution at nodes
    ue = uex(mesh.verts(:,1),mesh.verts(:,2));
    % assemble matrix and rhs for VEM and SSE-VEM
    [AVEM,ASSEVEM,b,out,A2] = SSEVEM(mesh,f);
    % SSE-VEM solution
    freeNode = setdiff(1:size(mesh.verts,1),mesh.bdNodes)';
    u = zeros(size(mesh.verts,1),1);
    u(freeNode) = ASSEVEM(freeNode,freeNode)\b(freeNode);
    % compute L2-error
    errSSEVEMhexa(pt+1) = getL2Error(out,uex,u); %(u-ue)'*A2*(u-ue);%
    % VEM solution
    u = zeros(size(mesh.verts,1),1);
    u(freeNode) = AVEM(freeNode,freeNode)\b(freeNode);
    % compute L2-error
    errVEMhexa(pt+1) = getL2Error(out,uex,u); % (u-ue)'*A2*(u-ue);%
    %  compute diam
    hhSSEVEMhexa(pt+1) = mesh.hmax;
    hhVEMhexa(pt+1) = mesh.hmax;
    pt = pt + 1;
end

%%
close all
figure
loglog(hhVEMvoro,errVEMvoro,'k.-'), hold on
loglog(hhVEMvoro,errSSEVEMvoro,'b.-')
loglog(hhVEMvoro,.05*hhVEMvoro.^2,'k.--'), hold off
grid on
set(gca,'FontSize',18)
legend({'VEM','SSE-VEM','$h^2$'},'Location','southeast','Interpreter','latex')
%axis([1e-2 .5 1e-5 1e-1])
%yticks([1e-5 1e-3 1e-1])

figure
loglog(hhVEMhexa,errVEMhexa,'k.-'), hold on
loglog(hhVEMhexa,errSSEVEMhexa,'b.-')
%loglog(hhVEM1,errLIVEM1,'r.-')
loglog(hhVEMhexa,.05*hhVEMhexa.^3,'k.--'), hold off
grid on
set(gca,'FontSize',18)
legend({'VEM','SSE-VEM','$h^2$'},'Location','southeast','Interpreter','latex')
%axis([1e-2 .5 1e-5 1e-1])
%yticks([1e-5 1e-3 1e-1])