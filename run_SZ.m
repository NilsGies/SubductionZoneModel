function run_SZ(sznumin)
%% clean up
clear -except all_data  close

if not(exist('sznumin','var'))
    sznumin=1:56;
end
%% defining your working directory
base_dir='/change/to/your/path/Final_Gies_et_al_2023'; 
in_folder='/SZ_Model';
output_directory=[base_dir '/data_output_ext_int1km'];

%% log calculation time
time.time_calc(1)= datetime('now');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% water saturation settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for water_sat=0 % water_sat=0 -> 2 % water in mantle ; water_sat=1 -> water saturated
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% selection of SZs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for wanted_SZ = sznumin %1:56   % 1) 1:8 2) 9:16 3) 17:24  4) 25:32 5) 33:40 6)  41:48 7)  49:56
        time.time_SZ(1,wanted_SZ) = datetime('now');

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% models with the following settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        % 1     hp02     Koma   (A-phaseK)   (mgsur)    BM
        % 2     hp02     PW     (A-phase)    (mgsur)    BM
        % 3     hp11     koma   (A-phAK11)   (mgsur11)  Tait
        % 4     hp11     PW     (A-phase)    (mgsur11)  Tait
        % 5     hp02     Koma   (A-phaseK)              BM
        % 6     hp02     PW     (A-phase)               BM
        % 7     hp11     koma   (A-phAK11)              Tait
        % 8     hp11     PW     (A-phase)               Tait
        %

        for model=[2,6,1]%[1:8]
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manual input section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%  copyfolder
            newdir=[base_dir in_folder '_' num2str(wanted_SZ) '_' num2str(model) '_' num2str(water_sat)];
            if exist(newdir,"dir")
               disp(['Skipping Model:' num2str(model)])
                clearvars -except wanted_SZ water_sat matrix_interpolation crust_only model total_h2o_recycled_SZ_TG_Ma sz_names  output_directory time output_dir_H2O base_dir  in_folder newdir
                cd(base_dir)
                continue
            end
            copyfile([base_dir, in_folder],newdir)
            cd(newdir)

            time.time_model(1,model)= datetime('now');



            %in order to check the changing bulk rock composition in the input file
            write_copy=0;

            %auto selection of paths and points
            auto_select=1;

            %if auto select = 0, determine the number of paths and points to be
            %calculated
            manual_paths=2;
            manual_points=250;

            %how is the oxygen treated?
            %it can be internally controlled (both values = 0), with a buffer
            %(with_buffer=1, with_oxygen = 0) or oxygen can be constrained (with_buffer=0, with_oxygen = 1)

            %oxygen as constrained component?
            with_oxygen=0;
            %or with a buffer?
            with_buffer=0;

            %calculations with water fractionation?
            water_fractionation=1;

            %do you use interpolated matrices?
            matrix_interpolation=4;

            %solution models
            amphibole1='GlTrTsPg';
            amphibole2='Amph(DPW)';
            omphacite='Omph(HP)';
            garnet='Gt(HP)';
            feldspar='feldspar';
            chlorite='Chl(HP)';
            epidote='Ep(HP)';
            olivine='O(HP)';
            opx='Opx(HP)';
            phengite='Pheng(HP)';
            lawsonite='   ';
            zoisite='  ';
            water='H2O';
            antigorite='Atg(PN)';
            cloritoid='Ctd(HP)';
            solution1=' ';
            solution2=' ';

            %considered (non solution) phases to be stored (can be plotted later)
            add_phase1='stlb';
            add_phase2='gth';
            add_phase3='zo';
            add_phase4='prl';
            add_phase5='br';
            add_phase6='heu';
            add_phase7='law';
            add_phase8='pump';



            % now the phases that differ in the models

            if model == 1||model ==5
                %solution models
                aphase='A-phaseK';

                %considered (non solution) phases to be stored (can be plotted later)
                add_phase9='mgsur';
            end
            if model == 2||model ==6

                %solution models
                aphase='A-phase';

                %considered (non solution) phases to be stored (can be plotted later)
                add_phase9='mgsur';

            end

            if model == 3||model ==7

                %solution models
                aphase='A-phAK11';

                %considered (non solution) phases to be stored (can be plotted later)
                add_phase9='mgsur11';
            end

            if model == 4||model ==8

                %solution models
                aphase='A-phase';

                %considered (non solution) phases to be stored (can be plotted later)
                add_phase9='mgsur11';
            end

            phase_data_cell=get_phase_data_cell;
            phase_data_cell_read=[chlorite;amphibole2;amphibole1;garnet;feldspar;phengite;olivine;...
                omphacite;opx;lawsonite;epidote;water;antigorite;cloritoid;solution1;solution2;...
                phase_data_cell(:,2);aphase;add_phase9];

            %name of the vertex input file. This file is only needed for the file
            %structure and the excluded phases
            initial='single_orig2.dat';

            %excluded phases
            excl1='qfm';
            excl2='mthm';
            excl3='hen';
            excl4='frw';
            excl5='h2oL';
            excl6='acti';
            excl7='vsv';
            excl8='kals';
            excl9='wo';
            excl10=' ';
            excl11=' ';

            %name of the thermodynamic database

            % model==1
            database_setting{1}='Databases/hp02_koma_mgsur.dat';
            % model==2
            database_setting{2}='Databases/hp02_PW_mgsur.dat';
            % model==3
            database_setting{3}='Databases/hp11_koma_mgsur.dat';
            % model==4
            database_setting{4}='Databases/hp11_PW_mgsur.dat';
            % model==5
            database_setting{5}='Databases/hp02_koma.dat';
            % model==6
            database_setting{6}='Databases/hp02_PW.dat';
            % model==7
            database_setting{7}='Databases/hp11_koma.dat';
            % model==8
            database_setting{8}='Databases/hp11_PW.dat';


            %name of the solution model
            solution='solution_short.dat';


            %% now the bulk rock compositions for each layer

            % make the bulk rock matrices for each PT path H2O, CO2, SiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O
            % slab structure: 12xslab mantle, 5x gabbro, 3x sheeted dykes, 1x seds, 8x
            % hanging wall mantle
            slab.size_hanging_wall_mantle=8;    %number of nodes in the wedge mantle layer
            slab.size_sediments=1;                   %number of nodes in the sediment layer
            slab.size_MORB=8;                   %number of nodes in the MORB layer
            slab.size_dykes=3;
            slab.size_gabbros=5;
            slab.size_slab_mantle_1=6;          %number of nodes in the 1st (upper) slab mantle layer
            slab.size_slab_mantle_2=6;          %number of nodes in the 2nd slab mantle layer

            if water_sat==0 % wt% composition for 2% water in mantle
                %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O  O2
                bulk_slab_mantle_1=[2; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005
                %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_slab_mantle_2=[2; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005
                %            H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_gabbros=[1; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
                %          H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_dykes=[3; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
                %         H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_seds=[7; 0; 58.57; 11.91; 5.21; 2.48; 5.95; 2.43; 2.04; 0]; %GLOSS Plank and Langmuir 1998
                %                 H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_wedge_mantle=[0; 0; 44.9; 4.44; 8.03; 37.71; 3.54; 0.36; 0; 0]; %PUM Workman and Hart 2005
            end

            if water_sat==1 % wt% composition for water saturated mantle
                %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O  O2
                bulk_slab_mantle_1=[40; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005 dry
                %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_slab_mantle_2=[40; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005 wet
                %            H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_gabbros=[1; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
                %          H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_dykes=[3; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
                %         H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_seds=[7; 0; 58.57; 11.91; 5.21; 2.48; 5.95; 2.43; 2.04; 0]; %GLOSS Plank and Langmuir 1998
                %                 H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
                bulk_wedge_mantle=[0; 0; 44.9; 4.44; 8.03; 37.71; 3.54; 0.36; 0; 0]; %PUM Workman and Hart 2005
            end

            % defining the slab
            slab.hw_mtl_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB+slab.size_sediments+slab.size_hanging_wall_mantle;
            slab.hw_mtl_start=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB+slab.size_sediments+1;

            slab.seds_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB+slab.size_sediments;
            slab.seds_start=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB+1;

            slab.morb_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB;
            slab.morb_start=slab.size_slab_mantle_1+slab.size_slab_mantle_2+1;

            slab.dykes_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_MORB;
            slab.dykes_start=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_gabbros+1;

            slab.gabbros_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2+slab.size_gabbros;
            slab.gabbros_start=slab.size_slab_mantle_1+slab.size_slab_mantle_2+1;

            slab.mtl_2_end=slab.size_slab_mantle_1+slab.size_slab_mantle_2;
            slab.mtl_2_start=slab.size_slab_mantle_1+1;

            slab.mtl_1_end=slab.size_slab_mantle_1;
            slab.mtl_1_start=1;


            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end manual input section %%%%%%%%%%%%%%%%%%%%%%%%

            %make compositional matrix
            compo_matrix = [bulk_slab_mantle_1,bulk_slab_mantle_2,bulk_gabbros,bulk_dykes,bulk_seds,bulk_wedge_mantle]';

            %conversion wt% into mol
            %mol weights
            H2O_mol_weight=18.01528;
            CO2_mol_weight=44.0095;
            SiO2_mol_weight=60.0843;
            Al2O3_mol_weight=101.9613;
            FeO_mol_weight=71.8464;
            MgO_mol_weight=40.3044;
            CaO_mol_weight=56.0794;
            Na2O_mol_weight=61.9789;
            K2O_mol_weight=94.196;
            TiO2_mol_weight=79.866;
            O_mol_weight=15.9994;
            %

            mol_weight_vector = [H2O_mol_weight, CO2_mol_weight, SiO2_mol_weight, Al2O3_mol_weight, FeO_mol_weight, MgO_mol_weight, CaO_mol_weight, Na2O_mol_weight, K2O_mol_weight, O_mol_weight];
            compo_matrix_mol = compo_matrix./mol_weight_vector;

            %% checking for old files and cleaning up
            if isfile('single_mol.arf')
                delete('single_mol.arf');
            end

            if isfile('single_mol.blk')
                delete('single_mol.blk');
            end

            if isfile('single_mol.plt')
                delete('single_mol.plt');
            end

            if isfile('single_mol.prn')
                delete('single_mol.prn');
            end

            if isfile('single_mol.tof')
                delete('single_mol.tof');
            end

            if isfile('control_out_mol')
                delete('control_out_mol');
            end

            if isfile('fractionated_bulk.dat')
                delete('fractionated_bulk.dat');
            end

            if isfile('H2O.dat')
                delete('H2O.dat');
            end

            if isfile('PTpoint')
                delete('PTpoint');
            end

            delete('single_mol_*.txt');

            %make the names of the files to be opened
            SZ_num=num2str(wanted_SZ);
            P_matrix_name='slabpaths/overlap_matrix_P_';
            T_matrix_name='slabpaths/overlap_matrix_T_';
            X_matrix_name='slabpaths/overlap_matrix_X_';
            Y_matrix_name='slabpaths/overlap_matrix_Y_';
            D_matrix_name='slabpaths/overlap_matrix_D_';
            H_matrix_name='slabpaths/overlap_matrix_H_';

            if matrix_interpolation==1
                P_matrix_name='slabpaths/overlap_matrix_P_int_';
                T_matrix_name='slabpaths/overlap_matrix_T_int_';
                X_matrix_name='slabpaths/overlap_matrix_X_int_';
                Y_matrix_name='slabpaths/overlap_matrix_Y_int_';
                D_matrix_name='slabpaths/overlap_matrix_D_int_';
                H_matrix_name='slabpaths/overlap_matrix_H_int_';
            end

            if matrix_interpolation==2
                P_matrix_name='slabpaths/overlap_matrix_P_int2_';
                T_matrix_name='slabpaths/overlap_matrix_T_int2_';
                X_matrix_name='slabpaths/overlap_matrix_X_int2_';
                Y_matrix_name='slabpaths/overlap_matrix_Y_int2_';
                D_matrix_name='slabpaths/overlap_matrix_D_int2_';
                H_matrix_name='slabpaths/overlap_matrix_H_int2_';
            end

            if matrix_interpolation==3
                P_matrix_name='slabpaths/overlap_matrix_P_int3_';
                T_matrix_name='slabpaths/overlap_matrix_T_int3_';
                X_matrix_name='slabpaths/overlap_matrix_X_int3_';
                Y_matrix_name='slabpaths/overlap_matrix_Y_int3_';
                D_matrix_name='slabpaths/overlap_matrix_D_int3_';
                H_matrix_name='slabpaths/overlap_matrix_H_int3_';
            end


if matrix_interpolation==4
                P_matrix_name='slabpaths/overlap_matrix_P_ext4_';
                T_matrix_name='slabpaths/overlap_matrix_T_ext4_';
                X_matrix_name='slabpaths/overlap_matrix_X_ext4_';
                Y_matrix_name='slabpaths/overlap_matrix_Y_ext4_';
                D_matrix_name='slabpaths/overlap_matrix_D_ext4_';
                H_matrix_name='slabpaths/overlap_matrix_H_ext4_';
            end


            P_path_matrix_name=[P_matrix_name SZ_num];
            T_path_matrix_name=[T_matrix_name SZ_num];
            X_path_matrix_name=[X_matrix_name SZ_num];
            Y_path_matrix_name=[Y_matrix_name SZ_num];
            D_path_matrix_name=[D_matrix_name SZ_num];
            H_path_matrix_name=[H_matrix_name SZ_num];


            %load the P,T,X and Y matrices for the SZ
            try
                P_path_matrix1=dlmread(P_path_matrix_name);
                T_path_matrix1=dlmread(T_path_matrix_name);
                X_path_matrix1=dlmread(X_path_matrix_name);
                Y_path_matrix1=dlmread(Y_path_matrix_name);
                D_path_matrix1=dlmread(D_path_matrix_name);
                H_path_matrix1=dlmread(H_path_matrix_name);
            catch ME
                P_path_matrix1=dlmread([P_path_matrix_name,'.txt']);
                T_path_matrix1=dlmread([T_path_matrix_name,'.txt']);
                X_path_matrix1=dlmread([X_path_matrix_name,'.txt']);
                Y_path_matrix1=dlmread([Y_path_matrix_name,'.txt']);
                D_path_matrix1=dlmread([D_path_matrix_name,'.txt']);
                H_path_matrix1=dlmread([H_path_matrix_name,'.txt']);
            end
            % defining the input file name
            name_input_file='single_mol';

            %define the P-T path of the metamorphic evolution
            %name_path_pre='reverse_paths1_';

            %determine the size
            size_PTXY_matrix=size(P_path_matrix1);

            %now make a figure of the SZ structure

            %define the SZ parameters
            fid=fopen('SZ_list.txt');
            datasz=fread(fid);
            szs=char(datasz');
            fclose(fid);
            sz_numbers =  textscan(szs,'%s%s');
            sz_names= sz_numbers(2);
            sz_number=(sz_numbers{1}(wanted_SZ));
            sz_name=char(sz_names{:}(wanted_SZ));
            %define the size and put it into the figure title
            paths_amaout=num2str(size_PTXY_matrix(2));
            node_amaount=num2str(size_PTXY_matrix(1));
            paths_amount_pre='Number of paths = ';
            nodes_amount_pre=' Number of points = ';
            name_pre=' SZ name = ';
            number_pre=' SZ number = ';

            %title_name=[number_pre sz_number name_pre sz_name paths_amount_pre paths_amaout nodes_amount_pre node_amaount];


            % uncomment this if you want to display the structure of the chosen subduction zone
            % figure(666)
            % plot(X_path_matrix1(1:size(Y_path_matrix1,1),:),Y_path_matrix1(1:size(Y_path_matrix1,1),:),'bo','MarkerSize',3)
            %
            %
            % title(title_name)
            % xlabel('Distance from trench (km)')
            % ylabel('Depth (km)')
            % daspect([1 1 1])
            %


            %%
            %now we have to make it one col larger as we need a col on top of the model
            %to extract the released water
            P_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            T_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            X_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            Y_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            D_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            H_path_matrix2=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);
            %matrix for the relative volumes for control
            rel_vol_matrix=zeros(size_PTXY_matrix(1),size_PTXY_matrix(2)+1);

            %fill the zeros matrix
            P_path_matrix2(:,1:end-1)=P_path_matrix1(:,:);
            T_path_matrix2(:,1:end-1)=T_path_matrix1(:,:);
            X_path_matrix2(:,1:end-1)=X_path_matrix1(:,:);
            Y_path_matrix2(:,1:end-1)=Y_path_matrix1(:,:);
            D_path_matrix2(:,1:end-1)=D_path_matrix1(:,:);
            H_path_matrix2(:,1:end-1)=H_path_matrix1(:,:);

            P_path_matrix2(:,end)=P_path_matrix1(:,end);
            T_path_matrix2(:,end)=T_path_matrix1(:,end);
            X_path_matrix2(:,end)=X_path_matrix1(:,end);
            Y_path_matrix2(:,end)=Y_path_matrix1(:,end);
            D_path_matrix2(:,end)=D_path_matrix1(:,end);
            H_path_matrix2(:,end)=H_path_matrix1(:,end);

            P_path_matrix=P_path_matrix2;
            T_path_matrix=T_path_matrix2;
            X_path_matrix=X_path_matrix2;
            Y_path_matrix=Y_path_matrix2;
            D_path_matrix=D_path_matrix2;
            H_path_matrix=H_path_matrix2;

            %define the number of paths, i.e. the layers and the number of points in
            %each layer
            if auto_select==1
                number_paths=size_PTXY_matrix(2);  %number of available PT paths
                number_points=size_PTXY_matrix(1); %number of PT points in these paths
            else
                %manual input
                number_paths=manual_paths;  %number of available PT paths
                number_points=manual_points; %number of PT points in these paths
            end

      %      number_points=3;

            %the suffix must be kept as vertex insists on the .dat extension
            input_suff='.dat';

            % make compositional zeros matrices
            H=zeros(number_paths,1); %number of PT paths
            C=zeros(number_paths,1);
            Si=zeros(number_paths,1);
            Al=zeros(number_paths,1);
            Fe=zeros(number_paths,1);
            Mg=zeros(number_paths,1);
            Ca=zeros(number_paths,1);
            Na=zeros(number_paths,1);
            K=zeros(number_paths,1);
            O=zeros(number_paths,1);

            for o = slab.mtl_1_start:1:slab.mtl_1_end
                %slab mantle 1
                H(o)=compo_matrix_mol(1,1);
                C(o)=compo_matrix_mol(1,2);
                Si(o)=compo_matrix_mol(1,3);
                Al(o)=compo_matrix_mol(1,4);
                Fe(o)=compo_matrix_mol(1,5);
                Mg(o)=compo_matrix_mol(1,6);
                Ca(o)=compo_matrix_mol(1,7);
                Na(o)=compo_matrix_mol(1,8);
                K(o)=compo_matrix_mol(1,9);
                O(o)=compo_matrix_mol(1,10);

            end

            for o = slab.mtl_2_start:1:slab.mtl_2_end
                %slab mantle 2
                H(o)=compo_matrix_mol(2,1);
                C(o)=compo_matrix_mol(2,2);
                Si(o)=compo_matrix_mol(2,3);
                Al(o)=compo_matrix_mol(2,4);
                Fe(o)=compo_matrix_mol(2,5);
                Mg(o)=compo_matrix_mol(2,6);
                Ca(o)=compo_matrix_mol(2,7);
                Na(o)=compo_matrix_mol(2,8);
                K(o)=compo_matrix_mol(2,9);
                O(o)=compo_matrix_mol(2,10);

            end

            for o = slab.gabbros_start:1:slab.gabbros_end
                %gabbros
                H(o)=compo_matrix_mol(3,1);
                C(o)=compo_matrix_mol(3,2);
                Si(o)=compo_matrix_mol(3,3);
                Al(o)=compo_matrix_mol(3,4);
                Fe(o)=compo_matrix_mol(3,5);
                Mg(o)=compo_matrix_mol(3,6);
                Ca(o)=compo_matrix_mol(3,7);
                Na(o)=compo_matrix_mol(3,8);
                K(o)=compo_matrix_mol(3,9);
                O(o)=compo_matrix_mol(3,10);
            end

            for o = slab.dykes_start:1:slab.dykes_end
                %dykes
                H(o)=compo_matrix_mol(4,1);
                C(o)=compo_matrix_mol(4,2);
                Si(o)=compo_matrix_mol(4,3);
                Al(o)=compo_matrix_mol(4,4);
                Fe(o)=compo_matrix_mol(4,5);
                Mg(o)=compo_matrix_mol(4,6);
                Ca(o)=compo_matrix_mol(4,7);
                Na(o)=compo_matrix_mol(4,8);
                K(o)=compo_matrix_mol(4,9);
                O(o)=compo_matrix_mol(4,10);
            end

            for o = slab.seds_start:1:slab.seds_end
                %sediments
                H(o)=compo_matrix_mol(5,1);
                C(o)=compo_matrix_mol(5,2);
                Si(o)=compo_matrix_mol(5,3);
                Al(o)=compo_matrix_mol(5,4);
                Fe(o)=compo_matrix_mol(5,5);
                Mg(o)=compo_matrix_mol(5,6);
                Ca(o)=compo_matrix_mol(5,7);
                Na(o)=compo_matrix_mol(5,8);
                K(o)=compo_matrix_mol(5,9);
                O(o)=compo_matrix_mol(5,10);
            end

            for o = slab.hw_mtl_start:1:slab.hw_mtl_end
                %wedge mantle
                H(o)=compo_matrix_mol(6,1);
                C(o)=compo_matrix_mol(6,2);
                Si(o)=compo_matrix_mol(6,3);
                Al(o)=compo_matrix_mol(6,4);
                Fe(o)=compo_matrix_mol(6,5);
                Mg(o)=compo_matrix_mol(6,6);
                Ca(o)=compo_matrix_mol(6,7);
                Na(o)=compo_matrix_mol(6,8);
                K(o)=compo_matrix_mol(6,9);
                O(o)=compo_matrix_mol(6,10);
            end

            %% make the name
            name_input=[name_input_file input_suff];

            %make control zeros matrices
            control_water_solids=zeros(number_paths,number_points);
            control_bulkH=zeros(number_paths,number_points);
            control_free_h2o=zeros(number_paths+1,number_points+1);
            control_add_h2o=zeros(number_paths+1,number_points+1);

            %make matrix for the free water that flows upward
            add_h2o=zeros(number_paths+1,number_points+1);
            add_CO2=zeros(number_paths+1,number_points+1);

            free_h2o_matrix=zeros(number_paths+1,number_points+1);
            free_CO2_matrix=zeros(number_paths+1,number_points+1);

            %make matrix for the water bound in solids
            wt_h2o_solids=zeros(number_paths+1,number_points+1);
            wt_CO2_solids=zeros(number_paths+1,number_points+1);

            %make a phase matrix
            phase_matrix=zeros(number_paths,number_points,33);
            phase_matrix_mol=zeros(number_paths,number_points,33);

            phasewt_2=zeros(number_paths,number_points,numel(phase_data_cell_read));
            
            if exist([output_directory '/SZ_results/' 'all_data_model_' num2str(model) '_H2O_2pct_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ)) '.mat'],'file')
                clearvars -except wanted_SZ water_sat matrix_interpolation crust_only model total_h2o_recycled_SZ_TG_Ma sz_names output_directory time output_dir_H2O base_dir  in_folder newdir
                cd(base_dir)
                rmdir(newdir,'s')
                continue
            end
            
            %% %%%%%%%%%%%% start calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %measure time needed for calculations
            startLoop = tic;
            time.loop_log=zeros(number_paths,1);
            time.iteration_log=zeros(number_paths,number_points);
            %count number of vertex calculations that failed finding the minimum
            count=0;
            endIteration=0;
            
            for path_to_calculate=1:number_paths     %number of available PT paths
                % state the database
                if model>4
                    database=database_setting{model};
                elseif model<=4 && path_to_calculate<=slab.mtl_2_end
                    database=database_setting{model};
                elseif model<=4 && path_to_calculate<slab.hw_mtl_start
                    database=database_setting{model+4};
                else
                    database=database_setting{model};
                end
                %at the first node of each P-T path make the initial bulk rock
                %composition from the above data

                %copy the original file for the file structure
                copyfile(initial,name_input);

                %make the PT point for the calculation

                %first read the P,T,X and Y vectors from the P,T,X and Y matrices
                P_path=P_path_matrix(:,path_to_calculate);
                T_path=T_path_matrix(:,path_to_calculate);
                X_path=X_path_matrix(:,path_to_calculate);
                Y_path=Y_path_matrix(:,path_to_calculate);
                D_path=D_path_matrix(:,path_to_calculate);
                H_path=H_path_matrix(:,path_to_calculate);

                %determine the P-T path size
                size_path=size(P_path);

                %now we have to make the initial input file as this is changing
                %according to the slab layering
                %open the input file
                fid = fopen(name_input,'r+t');
                %read the bulk rock composition
                bulk = textscan(fid, '%s%f%f','headerlines',29);
                fclose(fid);
                bulk_comp=bulk{3};

                H(path_to_calculate)=H(path_to_calculate)+add_h2o(path_to_calculate,1);
                C(path_to_calculate)=C(path_to_calculate)+add_CO2(path_to_calculate,1);

                %prepare the bulk for printing
                bulkH=sprintf('%s  %s  %s\n','H2O','1',num2str(H(path_to_calculate)));
                bulkC=sprintf('%s  %s  %s\n','CO2','1',num2str(C(path_to_calculate)));
                bulkSi=sprintf('%s  %s  %s\n','SIO2','1',num2str(Si(path_to_calculate)));
                bulkAl=sprintf('%s  %s  %s\n','AL2O3','1',num2str(Al(path_to_calculate)));
                bulkFe=sprintf('%s  %s  %s\n','FEO','1',num2str(Fe(path_to_calculate)));
                bulkMg=sprintf('%s  %s  %s\n','MGO','1',num2str(Mg(path_to_calculate)));
                bulkCa=sprintf('%s  %s  %s\n','CAO','1',num2str(Ca(path_to_calculate)));
                bulkNa=sprintf('%s  %s  %s\n','NA2O','1',num2str(Na(path_to_calculate)));
                bulkK=sprintf('%s  %s  %s\n','K2O','1',num2str(K(path_to_calculate)));
                if with_oxygen==1 %| with_buffer==1
                    bulkO=sprintf('%s  %s  %s\n','O2 ','1',num2str(O));
                else
                    bulkO=sprintf('%s  %s  %s\n','  ',' ','  ');
                end

                %now re-open the in_file
                fid = fopen(name_input,'r+t');
                %print new data into the input file
                fprintf(fid,'%s\r\n',database);
                fprintf(fid,'%s\r\n','print');
                fprintf(fid,'%s\r\n','plot');
                fprintf(fid,'%s\r\n',solution);
                fprintf(fid,'%s\r\n','plop');
                fprintf(fid,'%s\r\n','perplex_option.dat');
                fprintf(fid,'%s\r\n','10');
                fprintf(fid,'%s\r\n','PTpoint');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','15');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s\r\n','5');
                fprintf(fid,'%s\r\n','2');
                fprintf(fid,'%s\r\n','0');
                fprintf(fid,'%s %s %s %s %s\r\n','0.00000','0.00000','0.00000','0.00000','0.00000');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','begin','thermodynamic','component','list');
                fprintf(fid,'%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n',bulkH,bulkC,bulkSi,bulkAl,bulkFe,bulkMg,bulkCa,bulkNa,bulkK,bulkO);
                fprintf(fid,'%s %s %s %s\r\n','end','thermodynamic','component','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','begin','saturated','component','list');
                if with_buffer==1
                    fprintf(fid,'%s    %s  %s      %s      %s     %s\r\n','O2','0','0.00000','0.00000','0.00000','unconstrained amount');
                end
                fprintf(fid,'%s %s %s %s\r\n','end','saturated','component','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s %s\r\n','begin','saturated','phase','component','list');
                fprintf(fid,'%s %s %s %s %s\r\n','end','saturated','phase','component','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','begin','independent','potential/fugacity/activity','list');
                fprintf(fid,'%s %s %s %s\r\n','end','independent','potential/fugacity/activity','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','begin','excluded','phase','list');
                fprintf(fid,'%s\r\n',excl1);
                fprintf(fid,'%s\r\n',excl2);
                fprintf(fid,'%s\r\n',excl3);
                fprintf(fid,'%s\r\n',excl4);
                fprintf(fid,'%s\r\n',excl5);
                fprintf(fid,'%s\r\n',excl6);
                fprintf(fid,'%s\r\n',excl7);
                fprintf(fid,'%s\r\n',excl8);
                fprintf(fid,'%s\r\n',excl9);
                fprintf(fid,'%s\r\n',excl10);
                fprintf(fid,'%s\r\n',excl11);
                fprintf(fid,'%s %s %s %s\r\n','end','excluded','phase','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','begin','solution','phase','list');
                fprintf(fid,'%s\r\n',amphibole1);
                fprintf(fid,'%s\r\n',amphibole2);
                fprintf(fid,'%s\r\n',omphacite);
                fprintf(fid,'%s\r\n',garnet);
                fprintf(fid,'%s\r\n',feldspar);
                fprintf(fid,'%s\r\n',chlorite);
                fprintf(fid,'%s\r\n',epidote);
                fprintf(fid,'%s\r\n',opx);
                fprintf(fid,'%s\r\n',phengite);
                fprintf(fid,'%s\r\n',olivine);
                fprintf(fid,'%s\r\n',aphase);
                fprintf(fid,'%s\r\n',antigorite);
                fprintf(fid,'%s\r\n',cloritoid);
                fprintf(fid,'%s\r\n',solution1);
                fprintf(fid,'%s\r\n',solution2);
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%s %s %s %s\r\n','end','solution','phase','list');
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%f      %f     %f  %f      %f     %s %s %s %s %s %s\r\n',0,0,0,0,0,'max','p,','t,','xco2,','u1,','u2');
                fprintf(fid,'%f      %f     %f  %f      %f     %s %s %s %s %s %s\r\n',0,0,0,0,0,'max','p,','t,','xco2,','u1,','u2');
                fprintf(fid,'%f      %f     %f  %f      %f\r\n',0,0,0,0,0);
                fprintf(fid,'%s\r\n',' ');
                fprintf(fid,'%f  %f  %f  %f  %f\r\n',1,2,4,5,3);
                fclose(fid);

                % now the input with the correct bulk rock composition is ready for the calculations along the PT path

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin path loop %%%%%%%%%%%%%%%%%%%%%%%%%%
                for node_to_calculate=1:number_points   %number of nodes in each path
                    startIteration = tic;
                    nodes_done=path_to_calculate*number_points-number_points+node_to_calculate;
                    prozent_done=round(nodes_done/(number_paths*number_points)*100,2);


                    clc
                    if not(exist('calc_val','var')) || node_to_calculate==number_points 
                        calc_val=9;
                    end
                    calc_val=calc_val+1;
                    if calc_val==10
                        try
                            system("clear");
                        catch ME
                        end
                        disp(['Calculation Started: ' char(time.time_calc(1))])
                        disp(['Calculation SZ Started: ' char(time.time_SZ(1,wanted_SZ))])
                        disp(['Calculation Model Started: ' char(time.time_model(1,model))])
                        if water_sat==0
                            disp(['Calculating Model  # ', num2str(model), ' - with water = 2%'])
                        end
                        if water_sat==1
                            disp(['Calculating  Model # ', num2str(model), ' - water saturated'])
                        end

                        disp(['Calculating Subductionzone ', num2str(SZ_num), ' of 56: ',char(sz_name)])
                        disp(['Path (',num2str(path_to_calculate),'/',paths_amaout,')  ',' Node (',num2str(node_to_calculate),'/',node_amaount,')'])
                        disp(['Total (',num2str(nodes_done),'/',num2str(number_paths*number_points),')'])
                        disp([num2str(prozent_done),'%'])
                        disp(['Time for one node: ',num2str(endIteration)])
                        calc_val=0;
                    end
                    %make the PT point for calculation
                    PT_vertex(1,1)=P_path(node_to_calculate,1);
                    PT_vertex(1,2)=T_path(node_to_calculate,1);
                    dlmwrite('PTpoint',PT_vertex,' ');

                    %now write the PT point to file in order to make a PT matrix
                    Tmatrix=zeros(number_paths,size_path(1));
                    Pmatrix=zeros(number_paths,size_path(1));
                    Tmatrix(path_to_calculate,node_to_calculate)=PT_vertex(1,2);
                    Pmatrix(path_to_calculate,node_to_calculate)=PT_vertex(1,1);

                    if write_copy==1
                        control_pre1=num2str(path_to_calculate);
                        control_pre2=num2str(node_to_calculate);
                        name_copy_mol='control_input_mol.dat';
                        name_control_mol=[control_pre1 control_pre2 name_copy_mol];
       copyfile(name_input,name_control_mol);
                        
                    end

                    %% now run vertex and werami
                    
                    % with output printed to Command Window
                    %system('./PVK_batch_vertex_werami_new.sh')

                    % with suppressed output
                    [status, cmdout] = system('./PVK_batch_vertex_werami_mol_new.sh >NUL');


                    %% now extract the results
                    %first the phase abundances
                    %determine the name of the werami output file (this will be deleted later
                    %to keep the same name)
                    input_file_pre=name_input_file;
                    input_file_suff='_1.txt';

                    name_output=[input_file_pre input_file_suff];

                    %% %%%%%%%%%%%%%%% now the output read %%%%%%%%%%%%%%%%%%%%%

                    %% Read the output-file from vertex

                    %first I need the data
                    fid=fopen(name_output,'r+');
                    fid_error=0;
                    if fid==-1
                        fid_error=1;
                        name_output_subst='single_subst.txt';
                        %if compo_mol==0
                        %   fid=fopen('control_out','r+');
                        %else
                        fid=fopen('control_out_mol','r+');
                        %end

                    end
                    data=fread(fid);
                    s=char(data');
                    fclose(fid);
                    %to read the phase assemblage and the chemistry from the output file, first
                    %isolate the phases block utilizing keywords at the biginning and the end
                    key_size = strfind(s,'Solid+Melt');

                    if with_oxygen==0 && with_buffer==0
                        key1=strfind(s,'K2O');
                    else
                        key1=strfind(s,'O2');
                    end

                    %if minimisation failed open the last output file (control_out) and read
                    %the data from this file

                    if numel(key1)==0 %this is equivalent to an empty file in that case the control_out serves as output file
                        %if compo_mol==0
                        %    fid=fopen('control_out','r+');
                        %else
                        fid=fopen('control_out_mol','r+');
                        %end

                        data=fread(fid);
                        s=char(data');
                        fclose(fid);
                        key_size = strfind(s,'Solid+Melt');

                        if with_oxygen==0 && with_buffer==0
                            key_subst=strfind(s,'K2O');
                        else
                            key_subst=strfind(s,'O2');
                        end
                        count=count+1;
                        key_phases=key_subst(1);
                    else
                        key_phases=key1(1);
                    end

                    %first we extract the densities.
                    key_phases2=strfind(s,'Molar');
                    key_chem=strfind(s,'wt %');
                    key_molar=strfind(s,'Density(kg/m3)'); %14 characters

                    %here we must insert a switch in case water (or any fluid) is stable,
                    %because then there is this nasty "Sytem - fluid" keyword in the molar
                    %section

                    if numel(key_size)==0 %this is in case if there is no fluid (Solid+Melt) columns missing
                        key_chem=key_chem(2);
                        %read the small (no System - fluid row) molar section and write water
                        %as well as system density to matrix
                        key_molar2=strfind(s,'Seismic');
                        molar_block = s(key_molar+14:key_molar2-1);
                        [phases_molar, n_g_molar, enthalpy, entropy, volume, cp, alpha, beta, cp_cv, dens] =  strread(molar_block,'%s%f%f%f%f%f%f%f%f%f');
                        pos_molar1 = strmatch('System',phases_molar);
                        pos_molar2 = strmatch('H2O',phases_molar);
                        pos_molar3 = strmatch('CO2',phases_molar);
                        system_density = dens(pos_molar1(1));
                        phase_matrix(path_to_calculate,node_to_calculate,29)=system_density;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,29)=system_density;
                        %if water is present
                        if numel(pos_molar2)>0
                            h2o_density = dens(pos_molar2(1));
                        else
                            h2o_density = 0;
                            phase_matrix(path_to_calculate,node_to_calculate,28)=h2o_density;
                            phase_matrix_mol(path_to_calculate,node_to_calculate,28)=h2o_density;
                        end
                        %if CO2 is present
                        if numel(pos_molar3)>0
                            CO2_density = dens(pos_molar3(1));
                        else
                            CO2_density = 0;
                            phase_matrix(path_to_calculate,node_to_calculate,27)=CO2_density;
                            phase_matrix_mol(path_to_calculate,node_to_calculate,27)=CO2_density;
                        end

                    else
                        key_chem=key_chem(3);
                        %now again the molar block, but in the last row there is System - fluid
                        %as three strings, thus we have to skip the block into phases plus
                        %system and system - fluid rows

                        key_molar2=strfind(s,'System');
                        key_molar2=key_molar2(2);
                        key_molar3=strfind(s,'Seismic');
                        %phases
                        molar_block1 = s(key_molar+14:key_molar2-1);
                        %system - fluid
                        molar_block2 = s(key_molar2:key_molar3-1);
                        %phases: first water density
                        [phases_molar1, n_g_molar1, enthalpy1, entropy1, volume1, cp1, alpha1, beta1, cp_cv1, dens1] =  strread(molar_block1,'%s%f%f%f%f%f%f%f%f%f');
                        pos_molar2 = strmatch('H2O',phases_molar1);
                        h2o_density = dens1(pos_molar2(1));
                        phase_matrix(path_to_calculate,node_to_calculate,28)=h2o_density;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,28)=h2o_density;
                        %now for CO2
                        pos_molar3 = strmatch('CO2',phases_molar1);
                        %is CO2 present?
                        if numel(pos_molar3)>0
                            CO2_density = dens1(pos_molar3(1));
                        else
                            CO2_density = 0;
                        end
                        phase_matrix(path_to_calculate,node_to_calculate,27)=CO2_density;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,27)=CO2_density;
                        %system - fluid
                        [phases_molar2, hyphen, fluid_str, n_g_molar2, enthalpy2, entropy2, volume2, cp2, alpha2, beta2, cp_cv2, dens2] =  strread(molar_block2,'%s%s%s%f%f%f%f%f%f%f%f%f');
                        system_density = dens2;
                        phase_matrix(path_to_calculate,node_to_calculate,29)=system_density;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,29)=system_density;
                    end

                    key_chem2=strfind(s,'Other');

                    %read the phases and the chemistry blocks
                    %as O2 has only two letters we must change key phases + 3 into key phases
                    %+2
                    if with_oxygen==0 && with_buffer==0
                        phases_block = s(key_phases+4:key_phases2-1);
                    else
                        phases_block = s(key_phases+2:key_phases2-1);
                    end
                    chem_block = s(key_chem+4:key_chem2-1);

                    if with_oxygen==0 && with_buffer==0
                        %read data from phases blocks. First the phases block. Because there is an
                        %overflow character in the H2O, SiO2 and MgO in the Atg line these elements
                        %have to be read as char and then transformed
                        [phasewt,wtwt,volwt,molwt,mol,phaseH,phaseC,phaseSi,phaseAl,phaseFe,phaseMg,phaseCa,phaseNa,phaseK]=strread(phases_block,'%s%f%f%f%f%s%f%s%f%f%s%f%f%f');
                        phaseH = str2double(phaseH);
                        phaseSi = str2double(phaseSi);
                        phaseMg = str2double(phaseMg);
                    else
                        [phasewt,wtwt,volwt,molwt,mol,phaseH,phaseC,phaseSi,phaseAl,phaseFe,phaseMg,phaseCa,phaseNa,phaseK,phaseO]=strread(phases_block,'%s%f%f%f%f%s%f%s%f%f%s%f%f%f%f');
                        phaseH = str2double(phaseH);
                        phaseSi = str2double(phaseSi);
                        phaseMg = str2double(phaseMg);
                    end

                    % determin the size of the header
                    size_phasewt=size(phasewt);
                    % store the stable assemblage in an array for post checking
                    phasewt1{path_to_calculate,node_to_calculate}=phasewt;
                    pos_phasewt_2=cell(1,numel(phase_data_cell_read));
                    for n=1:numel(phase_data_cell_read)
                        pos_phasewt_2{n}=find(strcmp(phase_data_cell_read{n},phasewt));
                        if not(isempty(pos_phasewt_2{n}))
                            phasewt_2(path_to_calculate,node_to_calculate,n)=sum(wtwt(pos_phasewt_2{n}));
                        else
                            phasewt_2(path_to_calculate,node_to_calculate,n)=0;
                        end
                    end


                    %read the phase assemblage and fill the matrices
                    pos_phasewt{1}=strmatch(chlorite,phasewt);
                    pos_phasewt{2}=strmatch(amphibole2,phasewt);
                    pos_phasewt{3}=strmatch(amphibole1,phasewt);
                    pos_phasewt{4}=strmatch(garnet,phasewt);
                    pos_phasewt{5}=strmatch(feldspar,phasewt); %was zoesite
                    pos_phasewt{6}=strmatch(phengite,phasewt);
                    pos_phasewt{7}=strmatch(olivine,phasewt);
                    pos_phasewt{8}=strmatch(omphacite,phasewt);
                    pos_phasewt{9}=strmatch(opx,phasewt);
                    pos_phasewt{10}=strmatch(lawsonite,phasewt);
                    pos_phasewt{11}=strmatch(epidote,phasewt);
                    pos_phasewt{12}=strmatch(water,phasewt);
                    pos_phasewt{13}=strmatch(antigorite,phasewt);
                    pos_phasewt{14}=strmatch(aphase,phasewt);
                    pos_phasewt{15}=strmatch(cloritoid,phasewt);
                    pos_phasewt{16}=strmatch(solution1,phasewt);
                    pos_phasewt{17}=strmatch(solution2,phasewt);
                    %additional phases
                    pos_phasewt{18}=strmatch(add_phase1,phasewt); %add_phase1='   ' %add_phase1='phA'
                    pos_phasewt{19}=strmatch(add_phase2,phasewt); %add_phase2='gth'
                    pos_phasewt{20}=strmatch(add_phase3,phasewt); %add_phase3='zo'
                    pos_phasewt{21}=strmatch(add_phase4,phasewt); %add_phase4='prl'
                    pos_phasewt{22}=strmatch(add_phase5,phasewt); %add_phase6='heu'
                    pos_phasewt{23}=strmatch(add_phase6,phasewt);
                    pos_phasewt{24}=strmatch(add_phase7,phasewt); %add_phase7='law'
                    pos_phasewt{25}=strmatch(add_phase8,phasewt);
                    pos_phasewt{26}=strmatch(add_phase9,phasewt);%mgsur

                    %now put the phase abundances into the phase_matrix

                    for k=1:26
                        %determine the size of each position (solvus phases)
                        size_position=size(pos_phasewt{k});
                        sp=size_position(1);
                        %make a container matrix for the result of each phase and the sum
                        result=zeros(1,sp);
                        result_mol=zeros(1,sp);
                        sum_result=zeros(1,26);
                        sum_result_mol=zeros(1,26);

                        %fill the container matrix with the wt% or mol of each position
                        for l = 1:sp
                            pos= pos_phasewt{k}(l);
                            result(1,l) = wtwt(pos);
                            result_mol(1,l) = mol(pos);
                            sum_result(k)=sum(result);
                            sum_result_mol(k)=sum(result_mol);
                        end

                        %write the phase assemblage to matrix
                        phase_matrix(path_to_calculate,node_to_calculate,k)=sum_result(k);
                        phase_matrix_mol(path_to_calculate,node_to_calculate,k)=sum_result_mol(k);

                    end

                    size_grt_pos=size(pos_phasewt{4});
                    size_h2o_pos=size(pos_phasewt{12});

                    % extract the water data for fractionation
                    if size_h2o_pos(1)>0
                        h2o_fract_matrix(path_to_calculate,node_to_calculate)=mol(pos_phasewt{12});
                        if water_fractionation==1
                            fractH=h2o_fract_matrix(path_to_calculate,node_to_calculate);
                        end
                    else
                        fractH=0;
                    end

                    %then the chemistry block. This block has either 3 or 5 cols depending if
                    %melt or fluid is stable
                    if numel(key_size)==0
                        [elem_name,mol_chem,molpc_chem,wt_chem]=strread(chem_block,'%s%f%f%f');
                    else
                        [elem_name,mol_chem,molpc_chem,wt_chem,molpc_chem_solid,wtpc_chem_solid]=strread(chem_block,'%s%f%f%f%f%f');
                    end

                    %determine the positions of H2O and CO2

                    pos_phases1 = strmatch('H2O',phasewt);
                    pos_phases2 = strmatch('CO2',phasewt);

                    pos_chem1 = strmatch('H2O',elem_name);
                    pos_chem2 = strmatch('CO2',elem_name);

                    %extract free fluid
                    free_h2o_wt=wtwt(pos_phases1);
                    free_h2o_vol=volwt(pos_phases1);
                    free_h2o_molpc=molwt(pos_phases1);
                    free_h2o_mol=mol(pos_phases1);

                    %for control
                    if free_h2o_wt>=0
                        control_free_h2o(path_to_calculate,node_to_calculate)=free_h2o_wt;
                    else
                        control_free_h2o(path_to_calculate,node_to_calculate)=0;
                    end

                    free_CO2_wt=wtwt(pos_phases2);
                    free_CO2_vol=volwt(pos_phases2);
                    free_CO2_molpc=molwt(pos_phases2);
                    free_CO2_mol=mol(pos_phases2);

                    %extract fluids in bulk rock chemistry
                    if numel(key_size)==0
                        bulk_h2o_mol=mol_chem(pos_chem1);
                        bulk_h2o_molpc=molpc_chem(pos_chem1);
                        bulk_h2o_wtpc=wt_chem(pos_chem1);
                        bulk_CO2_mol=mol_chem(pos_chem2);
                        bulk_CO2_molpc=molpc_chem(pos_chem2);
                        bulk_CO2_wtpc=wt_chem(pos_chem2);
                        phase_matrix_mol(path_to_calculate,node_to_calculate,30)=bulk_h2o_molpc;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,31)=bulk_CO2_molpc;
                    else
                        bulk_h2o_mol=0;
                        bulk_CO2_mol=0;
                        bulk_h2o_molpc=molpc_chem_solid(pos_chem1);
                        bulk_h2o_wtpc=wtpc_chem_solid(pos_chem1);
                        bulk_CO2_molpc=molpc_chem_solid(pos_chem2);
                        bulk_CO2_wtpc=wtpc_chem_solid(pos_chem2);
                        phase_matrix_mol(path_to_calculate,node_to_calculate,30)=bulk_h2o_molpc;
                        phase_matrix_mol(path_to_calculate,node_to_calculate,31)=bulk_CO2_molpc;
                    end

                    wt_h2o_solids(path_to_calculate,node_to_calculate)=bulk_h2o_wtpc;
                    wt_CO2_solids(path_to_calculate,node_to_calculate)=bulk_CO2_wtpc;
try
                    phase_matrix(path_to_calculate,node_to_calculate,32)=bulk_h2o_wtpc;
                    phase_matrix(path_to_calculate,node_to_calculate,33)=bulk_CO2_wtpc;
catch ME
0
end
                    %now read the rest of the bulk rock chemistry for the next calculation step
                    %water amount to be fractionated is the molar amount of water


                    pos_phases3 = strmatch('SIO2',elem_name);
                    pos_phases4 = strmatch('AL2O3',elem_name);
                    pos_phases5 = strmatch('FEO',elem_name);
                    pos_phases6 = strmatch('MGO',elem_name);
                    pos_phases7 = strmatch('CAO',elem_name);
                    pos_phases8 = strmatch('NA2O',elem_name);
                    pos_phases9 = strmatch('K2O',elem_name);
                    pos_phases10 = strmatch('O2',elem_name);

                    bulkH=mol_chem(pos_chem1)-fractH;

                    if bulkH<0
                        bulkH=0;
                    end
                    bulkC=mol_chem(pos_chem2);
                    bulkSi=mol_chem(pos_phases3);
                    bulkAl=mol_chem(pos_phases4);
                    bulkFe=mol_chem(pos_phases5);
                    bulkMg=mol_chem(pos_phases6);
                    bulkCa=mol_chem(pos_phases7);
                    bulkNa=mol_chem(pos_phases8);
                    bulkK=mol_chem(pos_phases9);
                    bulkO=mol_chem(pos_phases10);



                    %% %%%%%%%%%%%%%%%%% end of new output read %%%%%%%%%%%%%%%%%%%%%
                    % first calcs with wt%
                    if pos_phases1>0
                        wt_free_h2o=free_h2o_molpc;
                    else
                        wt_free_h2o=0;
                    end
                    %
                    if pos_phases2>0
                        wt_free_CO2=free_CO2_molpc;
                    else
                        wt_free_CO2=0;
                    end

                    % % here we determin where the water flows to (vertical)
                    X_PT_point=X_path_matrix(node_to_calculate,path_to_calculate);
                    next_X_vector=X_path_matrix(:,path_to_calculate+1);
                    differ=abs(next_X_vector(:)-X_PT_point);

                    %find closes X value
                    [~, idx]=min(differ);

                    % now we have to determine the relative volumes of the outflow and inflow nodes
                    % information of which is in the D_ and H_matrices
                    outflow_volume=D_path_matrix(node_to_calculate,path_to_calculate).*H_path_matrix(node_to_calculate,path_to_calculate);
                    inflow_volume=D_path_matrix(idx,path_to_calculate+1).*H_path_matrix(idx,path_to_calculate+1);

                    volume_ratio_flow=outflow_volume/inflow_volume;
                    %write this value to matrix for control
                    rel_vol_matrix(node_to_calculate,path_to_calculate)=volume_ratio_flow;


                    if free_h2o_wt>=0
                        add_h2o(path_to_calculate+1,idx)=add_h2o(path_to_calculate+1,idx)+(free_h2o_mol*volume_ratio_flow);
                        free_h2o_matrix(path_to_calculate,node_to_calculate)=free_h2o_wt;
                        free_h2o_matrix_mol(path_to_calculate,node_to_calculate)=free_h2o_mol;
                    else
                        add_h2o(path_to_calculate+1,idx)=add_h2o(path_to_calculate+1,idx)+0;
                        free_h2o_matrix(path_to_calculate,node_to_calculate)=0;
                        free_h2o_matrix_mol(path_to_calculate,node_to_calculate)=0;
                    end
                    if free_CO2_wt>=0
                        add_CO2(path_to_calculate+1,idx)=add_CO2(path_to_calculate+1,idx)+(free_CO2_mol*volume_ratio_flow);
                        free_CO2_matrix(path_to_calculate,node_to_calculate)=free_CO2_wt;
                        free_CO2_matrix_mol(path_to_calculate,node_to_calculate)=free_CO2_mol;
                    else
                        add_CO2(path_to_calculate+1,idx)=add_CO2(path_to_calculate+1,idx)+0;
                        free_CO2_matrix(path_to_calculate,node_to_calculate)=0;
                        free_CO2_matrix_mol(path_to_calculate,node_to_calculate)=0;
                    end
                    %end

                    %% now open the input file for manipulation, save the bulk rock composition and close it again
                    copyfile(initial,name_input);

                    bulkH2=bulkH+add_h2o(path_to_calculate,node_to_calculate+1);
                    bulkC2=bulkC+add_CO2(path_to_calculate,node_to_calculate+1);

                    control_water(path_to_calculate,node_to_calculate)=wt_free_h2o;
                    control_bulkH(path_to_calculate,node_to_calculate)=bulkH;
                    control_add(path_to_calculate,node_to_calculate)=add_h2o(path_to_calculate,node_to_calculate);
                    control_bulkH2(path_to_calculate,node_to_calculate)=bulkH2;

                    % %prepare the bulk for printing
                    % %with corrected wt%
                    bulkH=sprintf('%s  %s  %s\n','H2O','1',num2str(bulkH2));
                    bulkC=sprintf('%s  %s  %s\n','CO2','1',num2str(bulkC2));
                    bulkSi=sprintf('%s  %s  %s\n','SIO2','1',num2str(bulkSi));
                    bulkAl=sprintf('%s  %s  %s\n','AL2O3','1',num2str(bulkAl));
                    bulkFe=sprintf('%s  %s  %s\n','FEO','1',num2str(bulkFe));
                    bulkMg=sprintf('%s  %s  %s\n','MGO','1',num2str(bulkMg));
                    bulkCa=sprintf('%s  %s  %s\n','CAO','1',num2str(bulkCa));
                    bulkNa=sprintf('%s  %s  %s\n','NA2O','1',num2str(bulkNa));
                    bulkK=sprintf('%s  %s  %s\n','K2O','1',num2str(bulkK));
                    if with_oxygen==1 %| with_buffer==1
                        bulkO=sprintf('%s  %s  %s\n','O2 ','1',num2str(bulkO,4));
                    else
                        bulkO=sprintf('%s  %s  %s\n','   ',' ','  ');
                    end

                    %%

                    %now re-open the in_file
                    fid = fopen(name_input,'r+');
                    % print new data into the input file
                    fprintf(fid,'%s\r\n',database);
                    fprintf(fid,'%s\r\n','print');
                    fprintf(fid,'%s\r\n','plot');
                    fprintf(fid,'%s\r\n',solution);
                    fprintf(fid,'%s\r\n','plop');
                    fprintf(fid,'%s\r\n','perplex_option.dat');
                    fprintf(fid,'%s\r\n','10');
                    fprintf(fid,'%s\r\n','PTpoint');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','15');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s\r\n','5');
                    fprintf(fid,'%s\r\n','2');
                    fprintf(fid,'%s\r\n','0');
                    fprintf(fid,'%s %s %s %s %s\r\n','0.00000','0.00000','0.00000','0.00000','0.00000');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','begin','thermodynamic','component','list');
                    fprintf(fid,'%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n',bulkH,bulkC,bulkSi,bulkAl,bulkFe,bulkMg,bulkCa,bulkNa,bulkK,bulkO);
                    fprintf(fid,'%s %s %s %s\r\n','end','thermodynamic','component','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','begin','saturated','component','list');
                    if with_buffer==1
                        fprintf(fid,'%s    %s  %s      %s      %s     %s\r\n','O2','0','0.00000','0.00000','0.00000','unconstrained amount');
                    end
                    fprintf(fid,'%s %s %s %s\r\n','end','saturated','component','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s %s\r\n','begin','saturated','phase','component','list');
                    fprintf(fid,'%s %s %s %s %s\r\n','end','saturated','phase','component','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','begin','independent','potential/fugacity/activity','list');
                    fprintf(fid,'%s %s %s %s\r\n','end','independent','potential/fugacity/activity','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','begin','excluded','phase','list');
                    fprintf(fid,'%s\r\n',excl1);
                    fprintf(fid,'%s\r\n',excl2);
                    fprintf(fid,'%s\r\n',excl3);
                    fprintf(fid,'%s\r\n',excl4);
                    fprintf(fid,'%s\r\n',excl5);
                    fprintf(fid,'%s\r\n',excl6);
                    fprintf(fid,'%s\r\n',excl7);
                    fprintf(fid,'%s\r\n',excl8);
                    fprintf(fid,'%s\r\n',excl9);
                    fprintf(fid,'%s\r\n',excl10);
                    fprintf(fid,'%s\r\n',excl11);
                    fprintf(fid,'%s %s %s %s\r\n','end','excluded','phase','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','begin','solution','phase','list');
                    fprintf(fid,'%s\r\n',amphibole1);
                    fprintf(fid,'%s\r\n',amphibole2);
                    fprintf(fid,'%s\r\n',omphacite);
                    fprintf(fid,'%s\r\n',garnet);
                    fprintf(fid,'%s\r\n',feldspar);
                    fprintf(fid,'%s\r\n',chlorite);
                    fprintf(fid,'%s\r\n',epidote);
                    fprintf(fid,'%s\r\n',opx);
                    fprintf(fid,'%s\r\n',phengite);
                    fprintf(fid,'%s\r\n',olivine);
                    fprintf(fid,'%s\r\n',aphase);
                    fprintf(fid,'%s\r\n',antigorite);
                    fprintf(fid,'%s\r\n',cloritoid);
                    fprintf(fid,'%s\r\n',solution1);
                    fprintf(fid,'%s\r\n',solution2);
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%s %s %s %s\r\n','end','solution','phase','list');
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%f      %f     %f  %f      %f     %s %s %s %s %s %s\r\n',0,0,0,0,0,'max','p,','t,','xco2,','u1,','u2');
                    fprintf(fid,'%f      %f     %f  %f      %f     %s %s %s %s %s %s\r\n',0,0,0,0,0,'max','p,','t,','xco2,','u1,','u2');
                    fprintf(fid,'%f      %f     %f  %f      %f\r\n',0,0,0,0,0);
                    fprintf(fid,'%s\r\n',' ');
                    fprintf(fid,'%f  %f  %f  %f  %f\r\n',1,2,4,5,3);
                    fclose(fid);

                    if numel(key1)>0 && fid_error<1
                        copyfile(name_output,'control_out_mol');
                    end %this is the end for the numel(key1) if
                    delete(name_output);
                    fid_error=0;
                    endIteration = toc(startIteration);
                    time.iteration_log(path_to_calculate,node_to_calculate)=endIteration;
                end %this is the end for the PT path loop
            end %this is the end for the paths


            %% now make a list of all phases that are stable during the PT loop
            % the names of the phases are already stored in the phasewt1{path_to_calculate} array

            % make an empty phase vector to store the phase names. There must be a
            % string that can be found in order to find the position for the next entry
            % of additional phases

            %phases_all = sprintf('%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty','empty');
            empty = cellstr('empty');
            phase_names(1:30,1) = empty;

            %now we write the first phase assemblage into the phase_names matrix
            size_phase_wt1 = size(phasewt1{1});
            phase_names(1:size_phase_wt1(1)) = phasewt1{1};

            endIteration = toc(startIteration);
            %% %%%%%%%%%%% calculation end, now result extraction %%%%%%%%%%%%%%%%%%%
            %now write the phase matrices
            chl_matrix=phase_matrix(:,:,1);
            atg_matrix=phase_matrix(:,:,13);
            amph_matrix=phase_matrix(:,:,3)+phase_matrix(:,:,2);
            grt_matrix=phase_matrix(:,:,4);
            fsp_matrix=phase_matrix(:,:,5);
            br_matrix=phase_matrix(:,:,22);
            phng_matrix=phase_matrix(:,:,6);
            ol_matrix=phase_matrix(:,:,7);
            cpx_matrix=phase_matrix(:,:,8);
            opx_matrix=phase_matrix(:,:,9);
            law_matrix=phase_matrix(:,:,24);
            ep_matrix=phase_matrix(:,:,11);
            aphase_matrix=phase_matrix(:,:,14);
            ctd_matrix=phase_matrix(:,:,15);
            h2o_matrix=phase_matrix(:,:,12);
            zo_matrix=phase_matrix(:,:,20);
            mgsur_matrix=phase_matrix(:,:,26);


            names=matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(phase_data_cell_read));
            for n=1:numel(names)
                if sum(phasewt_2(:,:,n),'all')>0
                    stable_phases.(char(names{n}))=phasewt_2(:,:,n);
                end
            end


            endLoop = toc(startLoop);
            time.loop_log(path_to_calculate)=endLoop;
            %%


            log_h2o_matrix=zeros(number_paths,number_points);
            for path_to_calculate=1:number_paths

                for node_to_calculate=1:number_points
                    if free_h2o_matrix(path_to_calculate,node_to_calculate)==0
                        log_h2o_matrix(path_to_calculate,node_to_calculate)=-3;
                    else
                        log_h2o_matrix(path_to_calculate,node_to_calculate)=log10(free_h2o_matrix(path_to_calculate,node_to_calculate));
                    end
                end
            end

            %% calculate total water
            % data from readme_van_keken
            % defining the parameters
            size_sediments=slab.size_sediments;               %number of nodes in the sediment layer
            size_MORB=slab.size_MORB;                   %number of nodes in the MORB layer
            size_hanging_wall_mantle=slab.size_hanging_wall_mantle;    %number of nodes in the wedge mantle layer
            size_slab_mantle_1=slab.size_slab_mantle_1;          %number of nodes in the 1st (upper) slab mantle layer
            size_slab_mantle_2=slab.size_slab_mantle_2;

            load SZ_info

            hight_volcanics=0.6*1000;
            hight_dykes=1.4*1000;
            hight_gabbros=5*1000;
            hight_slab_mantle=11*1000;
            speed_m_ma=SZ_info.Speedmmyr(wanted_SZ)*(1000/10000000); % m/ma
            hight_MORB=(hight_volcanics+hight_dykes+hight_gabbros)/size_MORB;
            hight_slab_mantle_single_node=hight_slab_mantle/(size_slab_mantle_1+size_slab_mantle_1);
            hight_matrix=zeros(29,1);
            for i=size_sediments+1:size_MORB+size_sediments
                hight_matrix(i)=hight_MORB;
            end
            for i=size_sediments+size_MORB+1:size_sediments+size_MORB+size_hanging_wall_mantle+size_slab_mantle_1+size_slab_mantle_2+1
                hight_matrix(i)=hight_slab_mantle_single_node;
            end
            hight_matrix=(flipud(H_path_matrix1(end ,:)'));
            length_SZ_in_m=SZ_info.Lengthkm(wanted_SZ)*1000;
            hight_seds=SZ_info.Sedimentkm(wanted_SZ)*1000;
            hight_seds_single_node=hight_seds/size_sediments;
            for i=1:size_sediments
                hight_matrix(i)=hight_seds_single_node;
            end

            % calculating the volume of the cell
            % xspace


            min_matrix_X=min(min(X_path_matrix1));
            max_matrix_X=max(max(X_path_matrix1));


            points_in_matrix=size(X_path_matrix1,1);

            lin_matrix_X=linspace(min_matrix_X,max_matrix_X,points_in_matrix);
            X_dist_1=(lin_matrix_X(2) - lin_matrix_X(1));

            % calculating the flux and recycled material
            speed_mm_a=SZ_info.Speedmmyr(i);
            speed_m_a=speed_mm_a*0.001;

            system_density_matrix=phase_matrix(:,:,29)'; %(kg/m3)
            volume_cell=speed_m_a.*(H_path_matrix1(1:number_points,1:number_paths)*1000*length_SZ_in_m); % m^3/a

            mass_flux_sz_kg=(system_density_matrix.*volume_cell); %kg/a
            cell_mass_flux_sz_kg_recycled=mass_flux_sz_kg(end,:);%kg/a
            total_mass_flux_sz_kg_recycled=sum(cell_mass_flux_sz_kg_recycled); %kg/a
            wth2o_recycled=wt_h2o_solids(1:number_paths,number_points-1)'/100;


            %%
            free_h2o_flux_Tg_Ma= free_h2o_matrix(1:end-1,1:end-1)'/100.*mass_flux_sz_kg*1e-9*1000000; %Tg/Ma; %
            h2o_flux_Tg_Ma= wt_h2o_solids(1:end-1,1:end-1)'/100.*mass_flux_sz_kg*1e-9*1000000; %Tg/Ma; %
if isempty(h2o_flux_Tg_Ma)
0
end 
    %%
            h2o_flux_kg_a=wth2o_recycled.*cell_mass_flux_sz_kg_recycled; %kg/a
            h2o_flux_SZ_TG_a=h2o_flux_kg_a*1e-9; % Tg/a
            h2o_recycled_SZ_TG_a=h2o_flux_SZ_TG_a(end,:); % Tg/a

            h2o_flux_SZ_TG_Ma=h2o_flux_SZ_TG_a*1000000; %Tg/Ma
            
            h2o_recycled_SZ_TG_Ma=h2o_flux_SZ_TG_Ma(end,:);
            total_h2o_recycled_SZ_TG_Ma=sum(h2o_recycled_SZ_TG_Ma);

            h2o_flux_kg_Ma=h2o_flux_kg_a*1000000;
            total_h2o_flux_kg_Ma=sum(h2o_flux_kg_Ma);

            %% saving water data
            if water_sat==0
                all_water_data.(['SZ_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ)) '_model_' num2str(model) '_H2O_2pct']) = total_h2o_recycled_SZ_TG_Ma;
            end

            if water_sat==1
                all_water_data.(['SZ_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ)) '_model_' num2str(model) '_H2O_sat']) = total_h2o_recycled_SZ_TG_Ma;
            end

            time.time_model(2,model)= datetime('now');

            %% save all data in structure
            w=whos;
            w_check=strcmp({w.name}, 'all_data');
            for get_it_sortet=1:length(w)
                if w_check(get_it_sortet)==0

                    if water_sat==0
                        all_data.(char(['model_' num2str(model) '_H2O_2pct_' char(sz_names{1}(wanted_SZ)) '_' num2str(wanted_SZ)])).(w(get_it_sortet).name) = eval(w(get_it_sortet).name);

                    end

                    if water_sat==1
                        all_data.(char(['model_' num2str(model) '_H2O_sat_' char(sz_names{1}(wanted_SZ)) '_' num2str(wanted_SZ)])).(w(get_it_sortet).name) = eval(w(get_it_sortet).name);
                    end

                end
            end


            %% save mat file in case one loop fails
            output_dir_mat=[output_directory '/SZ_results/'];
            try
                if ~exist(output_dir_mat, 'dir')
                    mkdir(output_dir_mat)
                end
                for get_it_sortet=1:length(w)
                    if w_check(get_it_sortet)==0

                        data_out.(char(w(get_it_sortet).name)) = eval(w(get_it_sortet).name);
                    end
                end
                if water_sat==0
                    save(([output_dir_mat 'all_data_model_' num2str(model) '_H2O_2pct_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ))]), '-struct', 'data_out', '-v7.3')
                end

                if water_sat==1
                    save(([output_dir_mat 'all_data_model_' num2str(model) '_H2O_sat_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ))]), '-struct', 'data_out', '-v7.3')
                end
            catch ME 
                0
            end

            %%
            output_dir_H2O=[output_directory '/total_water_results_Tg_Ma/'];

            if ~exist(output_dir_H2O, 'dir')
                mkdir(output_dir_H2O)
            end

            if water_sat==0

                save (([output_dir_H2O 'SZ_' num2str(wanted_SZ.','%02d') '_'  char(sz_names{1}(wanted_SZ)) '_model_' num2str(model) '_H2O_2pct']), '-struct', 'all_water_data')
            end

            if water_sat==1
                save (([output_dir_H2O 'SZ_' num2str(wanted_SZ.','%02d') '_' char(sz_names{1}(wanted_SZ)) '_model_' num2str(model) '_H2O_sat']), '-struct', 'all_water_data')
                clear total_h2o_recycled_SZ_TG_Ma
            end

            %% clear data
            clearvars -except wanted_SZ water_sat matrix_interpolation crust_only model total_h2o_recycled_SZ_TG_Ma sz_names output_directory time output_dir_H2O base_dir  in_folder newdir
            cd(base_dir)

            rmdir(newdir,'s')

        end %this is the end for the Model loop
    end %this is the end for the SZ loop
end %this is the end for the Water settings loop
