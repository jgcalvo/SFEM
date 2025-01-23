%% example 5.1
% comparisson for FEM, SFEM and VEM applied to structured and unstructured
% meshes
% exact solution and rhs
uex = @(x,y) sin(pi*x).*sin(pi*y);
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);

%% FEM
% read meshes and number of files
files     = dir('./meshes/stru_tria*');
numFiles  = size(files,1); 
errLIFEM1 = nan(numFiles,1);
hhFEM1    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [AFEM,b]   = FEM(mesh,f);
    % FEM solution
    uFEM = zeros(size(mesh.verts,1),1);
    uFEM(inNode) = AFEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    %errLIFEM1(pt+1) = norm(uFEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');%getL2Error_tria(uex,mesh,uFEM);
    %uElem = reshape(uFEM,3,mesh.nElems); uElem = uElem';
    errLIFEM1(pt+1) = getL2errorFEM(mesh,uex,uFEM);
    % compute diam
    hhFEM1(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

% read meshes and number of files
files     = dir('./meshes/unst_tria*');
numFiles  = size(files,1); 
errLIFEM2 = nan(numFiles,1);
hhFEM2    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [AFEM,b]   = FEM(mesh,f);
    % FEM solution
    uFEM = zeros(size(mesh.verts,1),1);
    uFEM(inNode) = AFEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    %errLIFEM2(pt+1) = norm(uFEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');%getL2Error_tria(uex,mesh,uFEM);
    %uElem = reshape(uFEM,3,mesh.nElems); uElem = uElem';
    errLIFEM2(pt+1) = getL2errorFEM(mesh,uex,uFEM);    
    % compute diam
    hhFEM2(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

%% ES-FEM
% read meshes and number of files
files     = dir('./meshes/stru_tria*');
numFiles  = size(files,1); 
errLIESFEM1 = nan(numFiles,1);
hhESFEM1    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [ASFEM,b]   = ESFEM(mesh,f);
    % FEM solution
    uSFEM = zeros(size(mesh.verts,1),1);
    uSFEM(inNode) = ASFEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    [~,~,out]   = VEM(mesh,f);
    errLIESFEM1(pt+1) = getL2Error(out,uex,uSFEM);
    %norm(uSFEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');%getL2Error_tria(uex,mesh,uFEM);
    % compute diam
    hhESFEM1(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

% read meshes and number of files
files     = dir('./meshes/unst_tria*');
numFiles  = size(files,1); 
errLIESFEM2 = nan(numFiles,1);
hhESFEM2    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [ASFEM,b]   = ESFEM(mesh,f);
    % FEM solution
    uSFEM = zeros(size(mesh.verts,1),1);
    uSFEM(inNode) = ASFEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    [~,~,out]   = VEM(mesh,f);
    errLIESFEM2(pt+1) = getL2Error(out,uex,uSFEM);
    %norm(uSFEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');%getL2Error_tria(uex,mesh,uFEM);
    % compute diam
    hhESFEM2(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

%% VEM
% read meshes and number of files
files     = dir('./meshes/stru_sfem*');
numFiles  = size(files,1); 
errLIVEM1 = nan(numFiles,1);
hhVEM1    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [AVEM,b,out]   = VEM(mesh,f);
    % VEM solution
    uVEM = zeros(size(mesh.verts,1),1);
    uVEM(inNode) = AVEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    %errLIVEM1(pt+1) = norm(uVEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');
    errLIVEM1(pt+1) = getL2Error(out,uex,uVEM);
    % compute diam
    hhVEM1(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

% read meshes and number of files
files     = dir('./meshes/unst_sfem*');
numFiles  = size(files,1); 
errLIVEM2 = nan(numFiles,1);
hhVEM2    = nan(numFiles,1);
pt = 0;
for idFile = 1:numFiles
    % load mesh
    mesh = load([files(idFile).folder '/' files(idFile).name]); 
    mesh = mesh.mesh;
    % compute boundary nodes
    inNode = setdiff(1:size(mesh.verts,1),mesh.bdNode);
    % assemble matrix and rhs for VEM and solve for interior nodes
    [AVEM,b,out]   = VEM(mesh,f);
    % VEM solution
    uVEM = zeros(size(mesh.verts,1),1);
    uVEM(inNode) = AVEM(inNode,inNode)\b(inNode);
    % compute L2-errors
    %errLIVEM2(pt+1) = norm(uVEM-uex(mesh.verts(:,1),mesh.verts(:,2)),'inf');%getL2Error_tria(uex,mesh,uFEM);
    errLIVEM2(pt+1) = getL2Error(out,uex,uVEM);
    % compute diam
    hhVEM2(pt+1) = max(mesh.diameter);
    pt = pt + 1;
end

%% 
close all
figure
loglog(hhFEM1,errLIFEM1,'k.-'), hold on
loglog(hhESFEM1,errLIESFEM1,'b.-')
loglog(hhVEM1,errLIVEM1,'r.-')
loglog(hhFEM1,.1*hhFEM1.^2,'k.--'), hold off
grid on
set(gca,'FontSize',18)
legend({'FEM','ES-FEM','VEM','$h^2$'},'Location','southeast','Interpreter','latex')
axis([1e-2 .5 1e-5 1e-1])
yticks([1e-5 1e-3 1e-1])

figure
loglog(hhFEM2,errLIFEM2,'k.-'), hold on
loglog(hhESFEM2,errLIESFEM2,'b.-')
loglog(hhVEM2,errLIVEM2,'r.-')
loglog(hhFEM2,.1*hhFEM2.^2,'k.--'), hold off
grid on
set(gca,'FontSize',18)
legend({'FEM','ES-FEM','VEM','$h^2$'},'Location','southeast','Interpreter','latex')
axis([1e-2 .5 1e-5 1e-1])
yticks([1e-5 1e-3 1e-1])


function L2err = getL2errorFEM(mesh,uex,uh)
node = mesh.verts;
elem = mesh.elems;
% compute L2err by adding L2 approximations of integral with one-point
% quadrature rule (evaluation at barycenter)
% approximated solution at barycenters
elem = reshape(cell2mat(elem),3,[])';
v1 = node(elem(:,1),:);
v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:);
xx = v1-v2; zz = v3-v1;
area = 0.5*(-xx(:,1).*zz(:,2) + xx(:,2).*zz(:,1));
% exact solution at barycenters
centroid = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
% compute quadrature rule
L2err = sqrt(sum(area.*(uex(centroid(:,1),centroid(:,2))-(uh(elem(:,1))+uh(elem(:,2))+uh(elem(:,3)))/3).^2));
end