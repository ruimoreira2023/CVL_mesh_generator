close all;clear all;clc

%% Data input %%

% Fixed and variable parameters
h1=2; % layer 1 thickness
h2=0.254; % layer 2 thickness
h3=2; % layer 3 thickness
L = 300; % length of the beam
Ns = [10]; % number of indentation
ds = [0.5]; % depths of indentation
alphas_req = [60]*pi/180; % angles of indentation
alphas_req_str = [60]; % angles to display on file name
fracs = [0.25]; % fractions to compute length of indentation

% Mesh size
na=6; % nof elements x-direction in region A
nb=2; % " region B
nc=6; % " region C
n3=4; % nof elements y-direction in layer 3
n2=2; % " layer 2
n1=3; % " layer 1
nz=10; % nof elements z-direction in a single element extrude
extrudeLength = 20; % length of the extrude command
no_delta=false; % if it does not consider (true) constant thickness of 
                % layer2 in indentation

%% Start generation of input files %%
for iNs = 1:length(Ns)
    
    mod = L/Ns(iNs)*0.5; % length of the module to then reflect
    cs = mod*fracs; % pitch of indentation
    
    for ids = 1:length(ds)
        
        d = ds(ids);
        
        for ics = 1:length(cs)
            
            c = cs(ics);
            
            for ialphas = 1:length(alphas_req)
                
                % Compute regions length
                a=mod - c - d*tan(alphas_req(ialphas));  % length region A
                b=mod - c - a;  % length region B
                
                % Compute angle so that h2 is kept constant or not
                beta=pi/2 - alphas_req(ialphas);
                alpha=(pi-beta)/2;
                
                % Find x distance (delta) to keep or not h2 constant
                if ~no_delta
                    delta=h2/tan(alpha);
                else
                    delta=0;
                end
                
                % Create geometry using extreme coordinates
                
                % camada 3 - Matriz 2x2
                % regiao A3
                XA3=[0 a+delta; 0 a+delta];
                YA3=[h1+h2+h3 h1+h2+h3; h1+h2 h1+h2];
                %regiao B3
                XB3=[a+delta a+delta+b; a+delta a+delta+b];
                YB3=[h1+h2+h3 h1+h2+h3; h1+h2 h1+h2-d];
                %regiao C3
                XC3=[a+delta+b a+b+c; a+delta+b a+b+c];
                YC3=[h1+h2+h3 h1+h2+h3; h1+h2-d h1+h2-d];
                
                % camada 2 - Matrix 2x2
                %regiao A2
                XA2=[0 a+delta; 0 a];
                YA2=[h1+h2 h1+h2; h1 h1];
                %regiao B2
                XB2=[a+delta a+delta+b; a a+b];
                YB2=[h1+h2 h1+h2-d; h1 h1-d];
                %regiao C2
                XC2=[a+delta+b a+b+c; a+b a+b+c];
                YC2=[h1+h2-d h1+h2-d; h1-d h1-d];
                
                %camada 1 - Matrix 2x2
                %regiao A1
                XA1=[0 a; 0 a];YA1=[h1 h1; 0 0];
                %regiao B1
                XB1=[a a+b; a a+b];
                YB1=[h1 h1-d; 0 0];
                %regiao C1
                XC1=[a+b a+b+c; a+b a+b+c];
                YC1=[h1-d h1-d; 0 0];
                
                % RA3x=[XA3(1,:); XA3(2,:); XA3(:,1)'; XA3(:,2)'];
                % RA3y=[YA3(1,:); YA3(2,:); YA3(:,1)'; YA3(:,2)'];
                % RB3x=[XB3(1,:); XB3(2,:); XB3(:,1)'; XB3(:,2)'];
                % RB3y=[YB3(1,:); YB3(2,:); YB3(:,1)'; YB3(:,2)'];
                % RC3x=[XC3(1,:); XC3(2,:); XC3(:,1)'; XC3(:,2)'];
                % RC3y=[YC3(1,:); YC3(2,:); YC3(:,1)'; YC3(:,2)'];
                %
                % RA2x=[XA2(1,:); XA2(2,:); XA2(:,1)'; XA2(:,2)'];
                % RA2y=[YA2(1,:); YA2(2,:); YA2(:,1)'; YA2(:,2)'];
                % RB2x=[XB2(1,:); XB2(2,:); XB2(:,1)'; XB2(:,2)'];
                % RB2y=[YB2(1,:); YB2(2,:); YB2(:,1)'; YB2(:,2)'];
                % RC2x=[XC2(1,:); XC2(2,:); XC2(:,1)'; XC2(:,2)'];
                % RC2y=[YC2(1,:); YC2(2,:); YC2(:,1)'; YC2(:,2)'];
                %
                % RA1x=[XA1(1,:); XA1(2,:); XA1(:,1)'; XA1(:,2)'];
                % RA1y=[YA1(1,:); YA1(2,:); YA1(:,1)'; YA1(:,2)'];
                % RB1x=[XB1(1,:); XB1(2,:); XB1(:,1)'; XB1(:,2)'];
                % RB1y=[YB1(1,:); YB1(2,:); YB1(:,1)'; YB1(:,2)'];
                % RC1x=[XC1(1,:); XC1(2,:); XC1(:,1)'; XC1(:,2)'];
                % RC1y=[YC1(1,:); YC1(2,:); YC1(:,1)'; YC1(:,2)'];
                %
                % figure(1); clf;
                % hl_a3=line(RA3x,RA3y);hl_b3=line(RB3x,RB3y);hl_c3=line(RC3x,RC3y);
                % hl_a2=line(RA2x,RA2y);hl_b2=line(RB2x,RB2y);hl_c2=line(RC2x,RC2y);
                % hl_a1=line(RA1x,RA1y);hl_b1=line(RB1x,RB1y);hl_c1=line(RC1x,RC1y);
                
                
                %USING PATCH
                p=[1 2 4 3];
                % Matrix 4x2 {line - x1,2,3,4 / column y1,2)}
                PA3= [XA3(p)',YA3(p)'] ;PB3= [XB3(p)',YB3(p)'] ;PC3= [XC3(p)',YC3(p)'] ;
                PA2= [XA2(p)',YA2(p)'] ;PB2= [XB2(p)',YB2(p)'] ;PC2= [XC2(p)',YC2(p)'] ;
                PA1= [XA1(p)',YA1(p)'] ;PB1= [XB1(p)',YB1(p)'] ;PC1= [XC1(p)',YC1(p)'] ;
                % Create geometry with patch
                figure(1); clf;
                hP_a3=patch(PA3(:,1),PA3(:,2),[0.6 0.6 0.6]);hP_b3=patch(PB3(:,1),PB3(:,2),[0.6 0.6 0.6]);
                hP_c3=patch(PC3(:,1),PC3(:,2),[0.6 0.6 0.6]);
                hP_a2=patch(PA2(:,1),PA2(:,2),[0, 100, 255] / 255);hP_b2=patch(PB2(:,1),PB2(:,2),[0, 100, 255] / 255);
                hP_c2=patch(PC2(:,1),PC2(:,2),[0, 100, 255] / 255);
                hP_a1=patch(PA1(:,1),PA1(:,2),[0.6 0.6 0.6]);hP_b1=patch(PB1(:,1),PB1(:,2),[0.6 0.6 0.6]);
                hP_c1=patch(PC1(:,1),PC1(:,2),[0.6 0.6 0.6]);
                
                % Write the name of each region/layer to understand better
                text(mean(XA3,"all"),mean(YA3,"all"),'A3');
                text(mean(XA2,"all"),mean(YA2,"all"),'A2');
                text(mean(XA1,"all"),mean(YA1,"all"),'A1');
                text(mean(XB3,"all"),mean(YB3,"all"),'B3');
                text(mean(XB2,"all"),mean(YB2,"all"),'B2');
                text(mean(XB1,"all"),mean(YB1,"all"),'B1');
                text(mean(XC3,"all"),mean(YC3,"all"),'C3');
                text(mean(XC2,"all"),mean(YC2,"all"),'C2');
                text(mean(XC1,"all"),mean(YC1,"all"),'C1');
                
                % New figure with same geometry, for meshing purposes
                figure(2); clf;
                hP_a3=patch(PA3(:,1),PA3(:,2),'w');hP_b3=patch(PB3(:,1),PB3(:,2),'w');
                hP_c3=patch(PC3(:,1),PC3(:,2),'w');
                hP_a2=patch(PA2(:,1),PA2(:,2),'w');hP_b2=patch(PB2(:,1),PB2(:,2),'w');
                hP_c2=patch(PC2(:,1),PC2(:,2),'w');
                hP_a1=patch(PA1(:,1),PA1(:,2),'w');hP_b1=patch(PB1(:,1),PB1(:,2),'w');
                hP_c1=patch(PC1(:,1),PC1(:,2),'w');
                
                %% ------- mapped nodes
                % Create mesh
                
                % region A3 map
                [XmA3,YmA3] = mapping(XA3,YA3,na,n3);
                [XmA3_reflected, YmA3_reflected] = reflectModule(XmA3,YmA3,mod);
                % region B3 map
                [XmB3,YmB3] = mapping(XB3,YB3,nb,n3);
                [XmB3_reflected, YmB3_reflected] = reflectModule(XmB3,YmB3,mod);
                % region C3 map
                [XmC3,YmC3] = mapping(XC3,YC3,nc,n3);
                [XmC3_reflected, YmC3_reflected] = reflectModule(XmC3,YmC3,mod);
                
                % region A2 map
                [XmA2,YmA2] = mapping(XA2,YA2,na,n2);
                [XmA2_reflected, YmA2_reflected] = reflectModule(XmA2,YmA2,mod);
                % region B2 map
                [XmB2,YmB2] = mapping(XB2,YB2,nb,n2);
                [XmB2_reflected, YmB2_reflected] = reflectModule(XmB2,YmB2,mod);
                % region C2 map
                [XmC2,YmC2] = mapping(XC2,YC2,nc,n2);
                [XmC2_reflected, YmC2_reflected] = reflectModule(XmC2,YmC2,mod);
                
                % region A1 map
                [XmA1,YmA1] = mapping(XA1,YA1,na,n1);
                [XmA1_reflected, YmA1_reflected] = reflectModule(XmA1,YmA1,mod);
                % region B1 map
                [XmB1,YmB1] = mapping(XB1,YB1,nb,n1);
                [XmB1_reflected, YmB1_reflected] = reflectModule(XmB1,YmB1,mod);
                % region C13 map
                [XmC1,YmC1] = mapping(XC1,YC1,nc,n1);
                [XmC1_reflected, YmC1_reflected] = reflectModule(XmC1,YmC1,mod);
                
                % z-direction coordinates for normal and reflected regions
                ZmA1 = extrudeModule(XmA1,extrudeLength,nz);ZmA1_reflected = extrudeModule(XmA1_reflected,extrudeLength,nz);
                ZmB1 = extrudeModule(XmB1,extrudeLength,nz);ZmB1_reflected = extrudeModule(XmB1_reflected,extrudeLength,nz);
                ZmC1 = extrudeModule(XmC1,extrudeLength,nz);ZmC1_reflected = extrudeModule(XmC1_reflected,extrudeLength,nz);
                
                ZmA2 = extrudeModule(XmA2,extrudeLength,nz);ZmA2_reflected = extrudeModule(XmA2_reflected,extrudeLength,nz);
                ZmB2 = extrudeModule(XmB2,extrudeLength,nz);ZmB2_reflected = extrudeModule(XmB2_reflected,extrudeLength,nz);
                ZmC2 = extrudeModule(XmC2,extrudeLength,nz);ZmC2_reflected = extrudeModule(XmC2_reflected,extrudeLength,nz);
                
                ZmA3 = extrudeModule(XmA3,extrudeLength,nz);ZmA3_reflected = extrudeModule(XmA3_reflected,extrudeLength,nz);
                ZmB3 = extrudeModule(XmB3,extrudeLength,nz);ZmB3_reflected = extrudeModule(XmB3_reflected,extrudeLength,nz);
                ZmC3 = extrudeModule(XmC3,extrudeLength,nz);ZmC3_reflected = extrudeModule(XmC3_reflected,extrudeLength,nz);
                
                % Draw mesh
                figure(2);
                hold on
                plot(XmA3,YmA3,'k');plot(XmA3',YmA3','k');plot(XmA3_reflected,YmA3_reflected,'k');plot(XmA3_reflected',YmA3_reflected','k');
                plot(XmB3,YmB3,'k');plot(XmB3',YmB3','k');plot(XmB3_reflected,YmB3_reflected,'k');plot(XmB3_reflected',YmB3_reflected','k');
                plot(XmC3,YmC3,'k');plot(XmC3',YmC3','k');plot(XmC3_reflected,YmC3_reflected,'k');plot(XmC3_reflected',YmC3_reflected','k');
                plot(XmA2,YmA2,'r');plot(XmA2',YmA2','r');plot(XmA2_reflected,YmA2_reflected,'k');plot(XmA2_reflected',YmA2_reflected','k');
                plot(XmB2,YmB2,'r');plot(XmB2',YmB2','r');plot(XmB2_reflected,YmB2_reflected,'k');plot(XmB2_reflected',YmB2_reflected','k');
                plot(XmC2,YmC2,'r');plot(XmC2',YmC2','r');plot(XmC2_reflected,YmC2_reflected,'k');plot(XmC2_reflected',YmC2_reflected','k');
                plot(XmA1,YmA1,'b');plot(XmA1',YmA1','b');plot(XmA1_reflected,YmA1_reflected,'k');plot(XmA1_reflected',YmA1_reflected','k');
                plot(XmB1,YmB1,'b');plot(XmB1',YmB1','b');plot(XmB1_reflected,YmB1_reflected,'k');plot(XmB1_reflected',YmB1_reflected','k');
                plot(XmC1,YmC1,'b');plot(XmC1',YmC1','b');plot(XmC1_reflected,YmC1_reflected,'k');plot(XmC1_reflected',YmC1_reflected','k');
                
                % nodal matrix
                
                % A1->A2->A3->B1->B2->B3->C1->C2->C3
                
                % Compute total number of nodes
                nnode=2*((na+1)*(n1+n2+n3+3)+(nb+1)*(n1+n2+n3+3)+(nc+1)*(n1+n2+n3+3)); % number of nodes counting the reflected ones
                nnodeExtrude = (nz+1)*nnode;
                
                % Each node will have [index, x, y, z]
                NODES=zeros(nnodeExtrude,4); % NODES=[inode xnode ynode znode];
                NODES(:,1)=1:nnodeExtrude;
                % First pattern of NODES coordinates -
                NODES(1:nnode,2)=[XmA1(:);XmA2(:);XmA3(:);XmB1(:);XmB2(:);XmB3(:);XmC1(:);XmC2(:);XmC3(:);...
                    XmA1_reflected(:);XmA2_reflected(:);XmA3_reflected(:);XmB1_reflected(:);XmB2_reflected(:);XmB3_reflected(:);XmC1_reflected(:);XmC2_reflected(:);XmC3_reflected(:)];
                NODES(1:nnode,3)=[YmA1(:);YmA2(:);YmA3(:);YmB1(:);YmB2(:);YmB3(:);YmC1(:);YmC2(:);YmC3(:);...
                    YmA1_reflected(:);YmA2_reflected(:);YmA3_reflected(:);YmB1_reflected(:);YmB2_reflected(:);YmB3_reflected(:);YmC1_reflected(:);YmC2_reflected(:);YmC3_reflected(:)];
                
                ZmA1 = coordToColumn(ZmA1,nz);ZmA2 = coordToColumn(ZmA2,nz);ZmA3 = coordToColumn(ZmA3,nz);
                ZmB1 = coordToColumn(ZmB1,nz);ZmB2 = coordToColumn(ZmB2,nz);ZmB3 = coordToColumn(ZmB3,nz);
                ZmC1 = coordToColumn(ZmC1,nz);ZmC2 = coordToColumn(ZmC2,nz);ZmC3 = coordToColumn(ZmC3,nz);
                ZmA1_reflected = coordToColumn(ZmA1_reflected,nz);ZmA2_reflected = coordToColumn(ZmA2_reflected,nz);ZmA3_reflected = coordToColumn(ZmA3_reflected,nz);
                ZmB1_reflected = coordToColumn(ZmB1_reflected,nz);ZmB2_reflected = coordToColumn(ZmB2_reflected,nz);ZmB3_reflected = coordToColumn(ZmB3_reflected,nz);
                ZmC1_reflected = coordToColumn(ZmC1_reflected,nz);ZmC2_reflected = coordToColumn(ZmC2_reflected,nz);ZmC3_reflected = coordToColumn(ZmC3_reflected,nz);
                
                NODES(1:nnode,4)=[ZmA1(:,1);ZmA2(:,1);ZmA3(:,1);ZmB1(:,1);ZmB2(:,1);ZmB3(:,1);ZmC1(:,1);ZmC2(:,1);ZmC3(:,1);...
                    ZmA1_reflected(:,1);ZmA2_reflected(:,1);ZmA3_reflected(:,1);ZmB1_reflected(:,1);ZmB2_reflected(:,1);ZmB3_reflected(:,1);ZmC1_reflected(:,1);ZmC2_reflected(:,1);ZmC3_reflected(:,1)];
                
                % Next pattern of NODES coordinates (actually didn't need it)
                M =(nnodeExtrude-nnode)/nnode;
                N = 1;
                NODES(nnode+1:nnodeExtrude,2)=repmat([XmA1(:);XmA2(:);XmA3(:);XmB1(:);XmB2(:);XmB3(:);XmC1(:);XmC2(:);XmC3(:);...
                    XmA1_reflected(:);XmA2_reflected(:);XmA3_reflected(:);XmB1_reflected(:);XmB2_reflected(:);XmB3_reflected(:);XmC1_reflected(:);XmC2_reflected(:);XmC3_reflected(:)], [M,N]);
                NODES(nnode+1:nnodeExtrude,3)=repmat([YmA1(:);YmA2(:);YmA3(:);YmB1(:);YmB2(:);YmB3(:);YmC1(:);YmC2(:);YmC3(:);...
                    YmA1_reflected(:);YmA2_reflected(:);YmA3_reflected(:);YmB1_reflected(:);YmB2_reflected(:);YmB3_reflected(:);YmC1_reflected(:);YmC2_reflected(:);YmC3_reflected(:)], [M,N]);
                
                
                
                x = 1;
                for izs = 1:nz
                    if x == 1
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 2
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 3
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 4
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 5
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 6
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 7
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 8
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 9
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    elseif x == 10
                        NODES(x*nnode+1:(x+1)*nnode,4)=[ZmA1(:,izs+1);ZmA2(:,izs+1);ZmA3(:,izs+1);ZmA1_reflected(:,izs+1);ZmA2_reflected(:,izs+1);ZmA3_reflected(:,izs+1);...
                            ZmB1(:,izs+1);ZmB2(:,izs+1);ZmB3(:,izs+1);ZmB1_reflected(:,izs+1);ZmB2_reflected(:,izs+1);ZmB3_reflected(:,izs+1);...
                            ZmC1(:,izs+1);ZmC2(:,izs+1);ZmC3(:,izs+1);ZmC1_reflected(:,izs+1);ZmC2_reflected(:,izs+1);ZmC3_reflected(:,izs+1)];
                    end
                    x = x+1;
                end
                
                % Copy nodes of entire module until length of the beam
                NODES=copyNodes(NODES,nnodeExtrude,Ns(iNs),mod);
                
                
                
                % element conect matrix
                % A1->A2->A3->B1->B2->B3->C1->C2->C3
                
                % Compute total number of elements
                nelem=2*(na+nb+nc)*(n1+n2+n3);
                nelemExtrude=nz*nelem;
                
                % Each element will have a quadrilateral form
                ELEM=zeros(nelemExtrude,7); %ELEM=[ielem node1 node2 node3 node4 Mat Prop];
                %ELEM=[ielem node1 node2 node3 node4 Mat Prop node5 node6 node7 node8 ];
                ELEM(:,1)=1:nelemExtrude;
                ELEM(:,6:7)=ones(nelemExtrude,2); %prop = 1 & mat = 1
                
                %generate elements for A1
                [ELEM,posE,posN]=genElem(ELEM,na,n1,0,0,1,4);
                %generate elements for A2
                [ELEM,posE,posN]=genElem(ELEM,na,n2,posE,posN,2,5);
                %generate elements for A3
                [ELEM,posE,posN]=genElem(ELEM,na,n3,posE,posN,1,6);
                %generate elements for B1
                [ELEM,posE,posN]=genElem(ELEM,nb,n1,posE,posN,1,4);
                %generate elements for B2
                [ELEM,posE,posN]=genElem(ELEM,nb,n2,posE,posN,2,5);
                %generate elements for B3
                [ELEM,posE,posN]=genElem(ELEM,nb,n3,posE,posN,1,6);
                %generate elements for C1
                [ELEM,posE,posN]=genElem(ELEM,nc,n1,posE,posN,1,4);
                %generate elements for C2
                [ELEM,posE,posN]=genElem(ELEM,nc,n2,posE,posN,2,5);
                %generate elements for C3
                [ELEM,posE,posN]=genElem(ELEM,nc,n3,posE,posN,1,6);
                
                % Generate reflected elements
                ELEM = reflectElem(ELEM,nelem,nnode);
                
                % Generate extruded elements
                [ELEM,numE,numN]=genElemExtrude(ELEM,na,n1,0,0,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,na,n2,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,na,n3,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nb,n1,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nb,n2,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nb,n3,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nc,n1,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nc,n2,numE,numN,nnode);
                [ELEM,numE,numN]=genElemExtrude(ELEM,nc,n3,numE,numN,nnode);
                
                % Generate reflected extruded elements
                ELEM = reflectElemExtruded(ELEM,nelem,nnode);
                
                % Generate all other extruded elements
                ELEM=genElemExtrudeMult(ELEM,nz,nelem,nnode);
                
                % Copy extruded module
                ELEM=copyElem(ELEM,nelemExtrude,nnodeExtrude,Ns(iNs));
                
                %                 % check coincident nodes
                %                 [NODES,ELEM]=CheckCoincidentNodes(NODES,ELEM,0.01);
                %
                %                 % renumber nodes
                %                 [NODES,ELEM]=RenumberNodes(NODES,ELEM,'+');
                
                %                 % Write indexes of nodes and elements to make it easier to understand
                %                 for inode=1:size(NODES,1)
                %                     text(NODES(inode,2),NODES(inode,3),NODES(inode,4),num2str(NODES(inode,1)),'fontsize',6);
                %                 end
                %                 for ielem=1:size(ELEM,1)
                %                     node_e=ELEM(ielem,2:5);[pnode_e,dum]=find(NODES(:,1)==node_e);
                %                     xloc=mean(NODES(pnode_e,2));
                %                     yloc=mean(NODES(pnode_e,3));
                %                     text(xloc, yloc, num2str(ielem),'fontsize',8,'Color','r')
                %                 end
                
                % Convert NODES matrix to write nastran format
                CP = zeros(size(NODES,1),1);
                NODES = [NODES(:,1) CP NODES(:,2)*10^(-3) NODES(:,3)*10^(-3) NODES(:,4)*10^(-3)];
                
                % Convert ELEM matrix to write nastran format
                %                 ELEM = [ELEM(1:2*nnode,1) ELEM(1:2*nnode,7) ELEM(1:2*nnode,2:5) ELEM(1:2*nnode,8:11)];
                ELEM = [ELEM(:,1) ELEM(:,7) ELEM(:,2:5) ELEM(:,8:11)];
                
                % Grid and CHEXA writing size
                [rNodes, cNodes] = size(NODES);
                [rElem, cElem] = size(ELEM);
                
                % Columns 2 to 7 writing of NODES
                col2 = [NODES(:,1)];
                col3 = [NODES(:,2)];
                col4 = [NODES(:,3)];
                col5 = [NODES(:,4)];
                col6 = [NODES(:,5)];
                col7 = [zeros(rNodes,1)];
                % Columns 2 to 11 writing of ELEM
                c2 = [ELEM(:,1)];
                c3 = [ELEM(:,2)];
                c4 = [ELEM(:,3)];
                c5 = [ELEM(:,4)];
                c6 = [ELEM(:,5)];
                c7 = [ELEM(:,6)];
                c8 = [ELEM(:,7)];
                c9 = [ELEM(:,8)];
                c10 = [ELEM(:,9)];
                c11 = [ELEM(:,10)];
                
                % Concatenate
                nastranFormatNODES = [col2,col3,col4,col5,col6,col7];
                nastranFormatELEM = [c2,c3,c4,c5,c6,c7,c8,c9,c10,c11];
                
                %%%%%% NASTRAN %%%%%%
                
                fnm = fullfile("newPreamble.dat");
                fidin = fopen(fnm,'r'); %opens the file for reading the Nastran Header format (Header_Nastran)
                folder = pwd;
%                               
                baseFileName=[num2str(Ns(iNs)), 'ind_d', num2str(d), '_c', num2str(c), '_ang', num2str(alphas_req_str(ialphas)), '_', num2str(h2), 'mm_toImport.dat'];
                fileName = fullfile(folder, baseFileName); %creates a file to import to Femap (Mesh_Matlab)
                fidout= fopen(fileName,'w'); %opens the file Mesh_Matlab for writing
                
                while true
                    thisline = fgets(fidin); %copy lines from Nastran Header
                    if ~ischar(thisline); break; end   %end of file copy
                    fwrite(fidout, thisline); %paste lines from Nastran Header into Mesh_Matlab
                end
                
                %Print new line to start format
                fprintf(fidout, '\n');
                
                %Print Grid
                % formatSpec1= 'GRID  \t %1.0f \t %1.0f \t %6.6f \t %6.6f \t %1.0f \t %1.0f\n';
                % formatSpec1= 'GRID  %10.0f \t %5.0f %12g %13g  %7.0f  %7.0f\n';
                formatSpec1= 'GRID,%1.0f,%1.0f,%g,%g,%g,%g\n';
                Mesh_grid=nastranFormatNODES(1:rNodes,:);
                fprintf(fidout,formatSpec1,Mesh_grid');
                
                %Print CHEXA
                % formatSpec2= 'CQUAD4 \t %1.0f \t %1.0f \t %1.6f \t %1.6f \t %1.0f \t %1.0f\n';
                % formatSpec2= 'CQUAD4  %8.0f \t %5.0f %12g %13g %8.0f  %7.0f\n';
                formatSpec2= 'CHEXA,%1.0f,%1.0f,%g,%g,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f\n';
                Mesh_chexa=nastranFormatELEM(1:end,:);
                fprintf(fidout,formatSpec2,Mesh_chexa');
                
                %Print enddata
                fprintf(fidout, '%s\n', 'ENDDATA'); %prints the last row code on Mesh_Matlab
                
                fclose(fidout); %close Mesh_Matlab
                fclose(fidin); %close Header_Nastran
                
                
            end
            
        end
        
    end
    
end
