results_data=dir('/change/to/your/path/Final_Gies_et_al_2023/data_output_ext_int1km/SZ_results/');
results_data=results_data(not([startsWith({results_data.name},'.')]) & not([results_data.isdir]));
 
load SZ_Model/SZ_info.mat
Data.Input.SZ_info=SZ_info;
Data.Input.SZ_names=SZ_info.Location;


for n2import=1:numel(results_data)
    
 nameSPLIT=split(results_data(n2import).name,'_');
SZ_num=str2double(nameSPLIT{7});
model=str2double(nameSPLIT{4});
    if isfield(Data,'Output') && isfield(Data.Output,'stable_phases') &&...
          numel( Data.Output)>=SZ_num &&...
        numel(Data.Output(SZ_num).stable_phases)>=model && not(isempty(Data.Output(SZ_num).stable_phases{model}))
    disp([num2str(n2import) ' of ' num2str(numel(results_data)) ' Skipped ' num2str(SZ_num) ' '  results_data(n2import).name ' model: ' num2str(model)] )
        continue
    end    
        load(fullfile(results_data(n2import).folder,results_data(n2import).name))

        SZ_num=str2double(SZ_num);

    disp([num2str(n2import) ' of ' num2str(numel(results_data)) ' Processing ' num2str(SZ_num) ' '  sz_name ' model: ' num2str(model)] )

    if model==1
        Data.Input.SZ_Variables(SZ_num).P_path_matrix1=P_path_matrix1;
        Data.Input.SZ_Variables(SZ_num).T_path_matrix1=T_path_matrix1;
        Data.Input.SZ_Variables(SZ_num).X_path_matrix1=X_path_matrix1;
        Data.Input.SZ_Variables(SZ_num).Y_path_matrix1=Y_path_matrix1;
        Data.Input.SZ_Variables(SZ_num).D_path_matrix1=D_path_matrix1;
        Data.Input.SZ_Variables(SZ_num).H_path_matrix1=H_path_matrix1;
    end


    %how is the oxygen treated?
    %it can be internally controlled (both values = 0), with a buffer
    %(with_buffer=1, with_oxygen = 0) or oxygen can be constrained (with_buffer=0, with_oxygen = 1)
    %oxygen as constrained component?
    Data.Input.Settings(model).with_oxygen=0;
    %or with a buffer?
    Data.Input.Settings(model).with_buffer=0;

    %calculations with water fractionation?
    Data.Input.Settings(model).water_fractionation=1;

  %  if SZ_num==56
        %solution models
        Data.Input.Settings(model).PerpleX.amphibole1='GlTrTsPg';
        Data.Input.Settings(model).PerpleX.amphibole2='Amph(DPW)';
        Data.Input.Settings(model).PerpleX.omphacite='Omph(HP)';
        Data.Input.Settings(model).PerpleX.garnet='Gt(HP)';
        Data.Input.Settings(model).PerpleX.feldspar='feldspar';
        Data.Input.Settings(model).PerpleX.chlorite='Chl(HP)';
        Data.Input.Settings(model).PerpleX.epidote='Ep(HP)';
        Data.Input.Settings(model).PerpleX.olivine='O(HP)';
        Data.Input.Settings(model).PerpleX.opx='Opx(HP)';
        Data.Input.Settings(model).PerpleX.phengite='Pheng(HP)';
        Data.Input.Settings(model).PerpleX.lawsonite='   ';
        Data.Input.Settings(model).PerpleX.zoisite='  ';
        Data.Input.Settings(model).PerpleX.water='H2O';
        Data.Input.Settings(model).PerpleX.antigorite='Atg(PN)';
        Data.Input.Settings(model).PerpleX.cloritoid='Ctd(HP)';
        Data.Input.Settings(model).PerpleX.solution1=' ';
        Data.Input.Settings(model).PerpleX.solution2=' ';

        %considered (non solution) phases to be stored (can be plotted later)
        Data.Input.Settings(model).PerpleX.add_phase1='stlb';
        Data.Input.Settings(model).PerpleX.add_phase2='gth';
        Data.Input.Settings(model).PerpleX.add_phase3='zo';
        Data.Input.Settings(model).PerpleX.add_phase4='prl';
        Data.Input.Settings(model).PerpleX.add_phase5='br';
        Data.Input.Settings(model).PerpleX.add_phase6='heu';
        Data.Input.Settings(model).PerpleX.add_phase7='law';
        Data.Input.Settings(model).PerpleX.add_phase8='pump';



        % now the phases that differ in the models
        if model == 1||model ==5
            %solution models
            Data.Input.Settings(model).PerpleX.aphase='A-phaseK';

            %considered (non solution) phases to be stored (can be plotted later)
            Data.Input.Settings(model).PerpleX.add_phase9='mgsur';
        end
        if model == 2||model ==6

            %solution models
            Data.Input.Settings(model).PerpleX.aphase='A-phase';

            %considered (non solution) phases to be stored (can be plotted later)
            Data.Input.Settings(model).PerpleX.add_phase9='mgsur';

        end

        if model == 3||model ==7

            %solution models
            Data.Input.Settings(model).PerpleX.aphase='A-phAK11';

            %considered (non solution) phases to be stored (can be plotted later)
            Data.Input.Settings(model).PerpleX.add_phase9='mgsur11';
        end

        if model == 4||model ==8

            %solution models
            Data.Input.Settings(model).PerpleX.aphase='A-phase';

            %considered (non solution) phases to be stored (can be plotted later)
            Data.Input.Settings(model).PerpleX.add_phase9='mgsur11';
        end




        %name of the vertex input file. This file is only needed for the file
        %structure and the excluded phases
        Data.Input.Settings(model).PerpleX.initial='single_orig2.dat';

        %excluded phases
        Data.Input.Settings(model).PerpleX.excl1='qfm';
        Data.Input.Settings(model).PerpleX.excl2='mthm';
        Data.Input.Settings(model).PerpleX.excl3='hen';
        Data.Input.Settings(model).PerpleX.excl4='frw';
        Data.Input.Settings(model).PerpleX.excl5='h2oL';
        Data.Input.Settings(model).PerpleX.excl6='acti';
        Data.Input.Settings(model).PerpleX.excl7='vsv';
        Data.Input.Settings(model).PerpleX.excl8='kals';
        Data.Input.Settings(model).PerpleX.excl9='wo';
        Data.Input.Settings(model).PerpleX.excl10=' ';
        Data.Input.Settings(model).PerpleX.excl11=' ';
        %name of the thermodynamic database

        % model==1
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp02_koma_mgsur.dat';
        % model==2
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp02_PW_mgsur.dat';
        % model==3
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp11_koma_mgsur.dat';
        % model==4
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp11_PW_mgsur.dat';
        % model==5
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp02_koma.dat';
        % model==6
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp02_PW.dat';
        % model==7
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp11_koma.dat';
        % model==8
        Data.Input.Settings(model).PerpleX.database_setting='Databases/hp11_PW.dat';


        %name of the solution model
        Data.Input.Settings(model).PerpleX.solution='solution_short.dat';
        water_sat=0;
        Data.Input.Global_SZ_Variables.slab.composition.Oxides={'H2O'  'CO2'  'SiO2'  'Al2O3'  'FeO'    'MgO'   'CaO'   'Na2O'  'K2O'  'O2'};

        if water_sat==0 % wt% composition for 2% water in mantle
            %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O  O2
            Data.Input.Global_SZ_Variables.slab.composition.bulk_slab_mantle_1=[2; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005
            %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_slab_mantle_2=[2; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005
            %            H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_gabbros=[1; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
            %          H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_dykes=[3; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
            %         H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_seds=[7; 0; 58.57; 11.91; 5.21; 2.48; 5.95; 2.43; 2.04; 0]; %GLOSS Plank and Langmuir 1998
            %                 H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_wedge_mantle=[0; 0; 44.9; 4.44; 8.03; 37.71; 3.54; 0.36; 0; 0]; %PUM Workman and Hart 2005
        end

        if water_sat==1 % wt% composition for water saturated mantle
            %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O  O2
            Data.Input.Global_SZ_Variables.slab.composition.bulk_slab_mantle_1=[40; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005 dry
            %                  H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_slab_mantle_2=[40; 0; 44.71; 3.98; 8.18; 38.73; 3.17; 0.13; 0; 0]; %DMM Workman and Hart 2005 wet
            %            H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_gabbros=[1; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
            %          H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_dykes=[3; 0; 49.51; 16.75; 8.05; 9.74; 12.5; 2.18; 0.065; 0]; %N-MORB Workman and Hart 2005
            %         H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_seds=[7; 0; 58.57; 11.91; 5.21; 2.48; 5.95; 2.43; 2.04; 0]; %GLOSS Plank and Langmuir 1998
            %                 H2O  CO2  SiO2  Al2O3  FeO    MgO   CaO   Na2O  K2O
            Data.Input.Global_SZ_Variables.slab.composition.bulk_wedge_mantle=[0; 0; 44.9; 4.44; 8.03; 37.71; 3.54; 0.36; 0; 0]; %PUM Workman and Hart 2005
        end
    %end

    %% now the bulk rock compositions for each layer

    % make the bulk rock matrices for each PT path H2O, CO2, SiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O
    % slab structure: 12xslab mantle, 5x gabbro, 3x sheeted dykes, 1x seds, 8x
    % hanging wall mantle
    %% SZ_num Global Inpuz Data

    Data.Input.Global_SZ_Variables.slab.size_hanging_wall_mantle=8;    %number of nodes in the wedge mantle layer
    Data.Input.Global_SZ_Variables.slab.size_sediments=1;                   %number of nodes in the sediment layer
    Data.Input.Global_SZ_Variables.slab.size_MORB=8;                   %number of nodes in the MORB layer
    Data.Input.Global_SZ_Variables.slab.size_dykes=3;
    Data.Input.Global_SZ_Variables.slab.size_gabbros=5;
    Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1=6;          %number of nodes in the 1st (upper) slab mantle layer
    Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2=6;          %number of nodes in the 2nd slab mantle layer



    % defining the slab
    Data.Input.Global_SZ_Variables.slab.hw_mtl_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB+Data.Input.Global_SZ_Variables.slab.size_sediments+Data.Input.Global_SZ_Variables.slab.size_hanging_wall_mantle;
    Data.Input.Global_SZ_Variables.slab.hw_mtl_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB+Data.Input.Global_SZ_Variables.slab.size_sediments+1;

    Data.Input.Global_SZ_Variables.slab.seds_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB+Data.Input.Global_SZ_Variables.slab.size_sediments;
    Data.Input.Global_SZ_Variables.slab.seds_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB+1;

    Data.Input.Global_SZ_Variables.slab.morb_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB;
    Data.Input.Global_SZ_Variables.slab.morb_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+1;

    Data.Input.Global_SZ_Variables.slab.dykes_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_MORB;
    Data.Input.Global_SZ_Variables.slab.dykes_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_gabbros+1;

    Data.Input.Global_SZ_Variables.slab.gabbros_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+Data.Input.Global_SZ_Variables.slab.size_gabbros;
    Data.Input.Global_SZ_Variables.slab.gabbros_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2+1;

    Data.Input.Global_SZ_Variables.slab.mtl_2_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+Data.Input.Global_SZ_Variables.slab.size_slab_mantle_2;
    Data.Input.Global_SZ_Variables.slab.mtl_2_start=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1+1;

    Data.Input.Global_SZ_Variables.slab.mtl_1_end=Data.Input.Global_SZ_Variables.slab.size_slab_mantle_1;
    Data.Input.Global_SZ_Variables.slab.mtl_1_start=1;

    %conversion wt% into mol
    %mol weights
    Data.Input.Global_SZ_Variables.H2O_mol_weight=18.01528;
    Data.Input.Global_SZ_Variables.CO2_mol_weight=44.0095;
    Data.Input.Global_SZ_Variables.SiO2_mol_weight=60.0843;
    Data.Input.Global_SZ_Variables.Al2O3_mol_weight=101.9613;
    Data.Input.Global_SZ_Variables.FeO_mol_weight=71.8464;
    Data.Input.Global_SZ_Variables.MgO_mol_weight=40.3044;
    Data.Input.Global_SZ_Variables.CaO_mol_weight=56.0794;
    Data.Input.Global_SZ_Variables.Na2O_mol_weight=61.9789;
    Data.Input.Global_SZ_Variables.K2O_mol_weight=94.196;
    Data.Input.Global_SZ_Variables.TiO2_mol_weight=79.866;
    Data.Input.Global_SZ_Variables.O_mol_weight=15.9994;
    %

    Data.Input.Global_SZ_Variables.mol_weight_vector = [Data.Input.Global_SZ_Variables.H2O_mol_weight, Data.Input.Global_SZ_Variables.CO2_mol_weight, Data.Input.Global_SZ_Variables.SiO2_mol_weight, Data.Input.Global_SZ_Variables.Al2O3_mol_weight, Data.Input.Global_SZ_Variables.FeO_mol_weight, Data.Input.Global_SZ_Variables.MgO_mol_weight, Data.Input.Global_SZ_Variables.CaO_mol_weight, Data.Input.Global_SZ_Variables.Na2O_mol_weight, Data.Input.Global_SZ_Variables.K2O_mol_weight, Data.Input.Global_SZ_Variables.O_mol_weight];
    Data.Input.Global_SZ_Variables.compo_matrix_mol = compo_matrix./mol_weight_vector;

    Data.Input.Global_SZ_Variables.phase_data_cell=get_phase_data_cell;
    Data.Input.Global_SZ_Variables.phase_data_cell_read=[Data.Input.Settings(model).PerpleX.chlorite;Data.Input.Settings(model).PerpleX.amphibole2;Data.Input.Settings(model).PerpleX.amphibole1;Data.Input.Settings(model).PerpleX.garnet;Data.Input.Settings(model).PerpleX.feldspar;Data.Input.Settings(model).PerpleX.phengite;Data.Input.Settings(model).PerpleX.olivine;...
        Data.Input.Settings(model).PerpleX.omphacite;Data.Input.Settings(model).PerpleX.opx;Data.Input.Settings(model).PerpleX.lawsonite;Data.Input.Settings(model).PerpleX.epidote;Data.Input.Settings(model).PerpleX.water;Data.Input.Settings(model).PerpleX.antigorite;Data.Input.Settings(model).PerpleX.cloritoid;Data.Input.Settings(model).PerpleX.solution1;Data.Input.Settings(model).PerpleX.solution2;...
        phase_data_cell(:,2);aphase;Data.Input.Settings(model).PerpleX.add_phase9];


%% results

Data.Output(SZ_num).stable_phases{model}=stable_phases;

%%Amph
Data.Output(SZ_num).VIP(model).amph_matrix=amph_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).amph_matrix=stable_phases.(char(matlab.lang.makeValidName(amphibole1)))+stable_phases.(char(matlab.lang.makeValidName(amphibole2)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).amph_matrix=[];
end

%%PhaseA
Data.Output(SZ_num).VIP(model).aphase_matrix=aphase_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).aphase_matrix=stable_phases.(char(matlab.lang.makeValidName(aphase)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).aphase_matrix=[];
end

%%Atg
Data.Output(SZ_num).VIP(model).atg_matrix=atg_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).atg_matrix=stable_phases.(char(matlab.lang.makeValidName(antigorite)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).atg_matrix=[];
end

%%Brucite
Data.Output(SZ_num).VIP(model).br_matrix=br_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).br_matrix=stable_phases.br;
catch
Data.Output(SZ_num).VIP_stable_phases(model).br_matrix=[];
end

%%Chlorite
Data.Output(SZ_num).VIP(model).chl_matrix=chl_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).chl_matrix=stable_phases.(char(matlab.lang.makeValidName(chlorite)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).chl_matrix=[];
end

%%Cloritoid
Data.Output(SZ_num).VIP(model).ctd_matrix=ctd_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).ctd_matrix=stable_phases.(char(matlab.lang.makeValidName(cloritoid)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).ctd_matrix=[];
end

%%CPX
Data.Output(SZ_num).VIP(model).cpx_matrix=cpx_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).cpx_matrix=stable_phases.(char(matlab.lang.makeValidName(omphacite)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).cpx_matrix=[];
end

%%Epidote
Data.Output(SZ_num).VIP(model).ep_matrix=ep_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).ep_matrix=stable_phases.(char(matlab.lang.makeValidName(epidote)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).ep_matrix=[];
end

%%Feldspar
Data.Output(SZ_num).VIP(model).fsp_matrix=fsp_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).fsp_matrix=stable_phases.(char(matlab.lang.makeValidName(feldspar)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).fsp_matrix=[];
end

%%grt
Data.Output(SZ_num).VIP(model).grt_matrix=grt_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).grt_matrix=stable_phases.(char(matlab.lang.makeValidName(grt_matrix)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).grt_matrix=[];
end

%%H2O
Data.Output(SZ_num).VIP(model).h2o_matrix=h2o_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).h2o_matrix=stable_phases.H2O;
catch
Data.Output(SZ_num).VIP_stable_phases(model).h2o_matrix=[];
end

%%law
Data.Output(SZ_num).VIP(model).law_matrix=law_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).law_matrix=stable_phases.law;
catch
Data.Output(SZ_num).VIP_stable_phases(model).law_matrix=[];
end

%%mgsur
Data.Output(SZ_num).VIP(model).mgsur_matrix=mgsur_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).mgsur_matrix=stable_phases.mgsur;
catch
Data.Output(SZ_num).VIP_stable_phases(model).mgsur_matrix=[];
end

%%olivine
Data.Output(SZ_num).VIP(model).ol_matrix=ol_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).ol_matrix=stable_phases.(char(matlab.lang.makeValidName(olivine)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).ol_matrix=[];
end

%%opx
Data.Output(SZ_num).VIP(model).opx_matrix=opx_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).opx_matrix=stable_phases.(char(matlab.lang.makeValidName(opx)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).opx_matrix=[];
end

%%phng
Data.Output(SZ_num).VIP(model).phng_matrix=phng_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).phng_matrix=stable_phases.(char(matlab.lang.makeValidName(phengite)));
catch
Data.Output(SZ_num).VIP_stable_phases(model).phng_matrix=[];
end
%%Zo
Data.Output(SZ_num).VIP(model).zo_matrix=zo_matrix;
try
Data.Output(SZ_num).VIP_stable_phases(model).zo_matrix=stable_phases.cz;
catch
Data.Output(SZ_num).VIP_stable_phases(model).zo_matrix=[];
end


try
Data.Output(SZ_num).VIP_stable_phases(model).ta_matrix=stable_phases.ta;
catch
Data.Output(SZ_num).VIP_stable_phases(model).ta_matrix=[];
end

try
Data.Output(SZ_num).system_density_matrix{model}=system_density_matrix;
Data.Output(SZ_num).mass_flux_sz_kg{model}=mass_flux_sz_kg;
Data.Output(SZ_num).volume_cell{model}=volume_cell;

Data.Output(SZ_num).free_h2o_flux_Tg_Ma{model}=free_h2o_flux_Tg_Ma;
Data.Output(SZ_num).free_h2o_matrix{model}=free_h2o_matrix;
Data.Output(SZ_num).free_h2o_matrix_mol{model}=free_h2o_matrix_mol;
Data.Output(SZ_num).h2o_flux_Tg_Ma{model}=h2o_flux_Tg_Ma;
Data.Output(SZ_num).h2o_fract_matrix{model}=h2o_fract_matrix;
Data.Output(SZ_num).h2o_matrix{model}=h2o_matrix;
Data.Output(SZ_num).wt_h2o_solids{model}=wt_h2o_solids;

catch  
0
end
end
disp('done')

for SZ_num=1:56
     for model=1:8
        if not(isfield(Data,'Output')) || numel(Data.Output)>SZ_num || model>size(Data.Output(SZ_num).system_density_matrix,2);
            continue
        end
       
        system_density_matrix=Data.Output(SZ_num).system_density_matrix{model};
        volume_cell=Data.Output(SZ_num).volume_cell{model};
        wt_h2o_solids=Data.Output(SZ_num).wt_h2o_solids{model};

        mass_flux_sz_kg=Data.Output(SZ_num).mass_flux_sz_kg{model};
        h2o_flux_Tg_Ma=Data.Output(SZ_num).h2o_flux_Tg_Ma{model};

        if isempty(h2o_flux_Tg_Ma)
    continue
        end


        h2o_flux_Tg_Ma_recycled{SZ_num}(:,model)=h2o_flux_Tg_Ma(end,:);

       %% water calculation check
        %% calculate total water
        % defining the parameters
        size_sediments=Data.Input.Global_SZ_Variables.slab.size_sediments;               %number of nodes in the sediment layer
       
        hight_matrix=(flipud(Data.Input.SZ_Variables(SZ_num).H_path_matrix1(end ,:)'))*1000; % m
        length_SZ_in_m=Data.Input.SZ_info.Lengthkm(SZ_num)*1000;  % m

        % Variable sediment hight for each SZ_num
        hight_seds=Data.Input.SZ_info.Sedimentkm(SZ_num)*1000; % m
        hight_seds_single_node=hight_seds/size_sediments;  % m

        hight_matrix(Data.Input.Global_SZ_Variables.slab.seds_start:Data.Input.Global_SZ_Variables.slab.seds_end)=hight_seds_single_node;



        % calculating the flux and recycled material
        speed_mm_a=Data.Input.SZ_info.Speedmmyr(SZ_num); % mm/a
        speed_m_a=speed_mm_a*0.001; % m/a

        % system_density_matrix=phase_matrix(:,:,29)'; %(kg/m3)
        [number_points,number_paths]=size(system_density_matrix);

        volume_cell=speed_m_a.*(Data.Input.SZ_Variables(SZ_num).H_path_matrix1(1:number_points,1:number_paths)*1000*length_SZ_in_m); % m^3/a
        mass_flux_sz_kg=(system_density_matrix.*volume_cell); %kg/a
        cell_mass_flux_sz_kg_recycled=mass_flux_sz_kg(end,:);%kg/a
        total_mass_flux_sz_kg_recycled=sum(cell_mass_flux_sz_kg_recycled); %kg/a
        wth2o_recycled=wt_h2o_solids(1:number_paths,number_points-1)'/100;


        %%
        h2o_flux_kg_a=  (wt_h2o_solids(1:end-1,1:end-1)'/100).*mass_flux_sz_kg; %kg/a; %
        h2o_flux_Tg_Ma= (wt_h2o_solids(1:end-1,1:end-1)'/100).*mass_flux_sz_kg*1e-9*1000000; %Tg/Ma; %

        free_h2o_matrix=Data.Output(SZ_num).free_h2o_matrix{model};
        free_h2o_flux_Tg_Ma= free_h2o_matrix(1:end-1,1:end-1)'/100.*mass_flux_sz_kg*1e-9*1000000; %Tg/Ma; %

        %%
        h2o_flux_kg_a_recycled=wth2o_recycled.*cell_mass_flux_sz_kg_recycled; %kg/a
        h2o_flux_kg_Ma=h2o_flux_kg_a_recycled*1000000;

        h2o_flux_SZ_TG_a_recycled=h2o_flux_kg_a_recycled*1e-9; % Tg/a
        h2o_flux_SZ_TG_a_recycled=h2o_flux_SZ_TG_a_recycled(end,:); % Tg/a

        h2o_flux_SZ_TG_Ma_recycled=h2o_flux_SZ_TG_a_recycled*1000000; %Tg/Ma

        total_h2o_recycled_SZ_TG_Ma=sum(h2o_flux_SZ_TG_Ma_recycled); %Tg/Ma

        Data.Output(SZ_num).global_h2o_flux_kg_Ma(model,1)=sum(h2o_flux_kg_Ma); %kg/Ma
        Data.Output(SZ_num).global_h2o_flux_Tg_Ma(model,1)=sum(total_h2o_recycled_SZ_TG_Ma); %Tg/Ma

        Data.Output(SZ_num).system_density_matrix{model}=system_density_matrix;
        Data.Output(SZ_num).mass_flux_sz_kg{model}=mass_flux_sz_kg;
        Data.Output(SZ_num).volume_cell{model}=volume_cell;

        Data.Output(SZ_num).free_h2o_flux_Tg_Ma{model}=free_h2o_flux_Tg_Ma;
        Data.Output(SZ_num).free_h2o_matrix{model}=free_h2o_matrix;

        Data.Output(SZ_num).h2o_flux_Tg_Ma{model}=h2o_flux_Tg_Ma;

        Data.Output(SZ_num).wt_h2o_solids{model}=wt_h2o_solids;
    end
end
%%
try

thx=[Data.Output(:).global_h2o_flux_Tg_Ma];
thx=[thx sum(thx,2)];
t=array2table(thx([1,2,6],:)','VariableNames',{'Dataset-1_with_Mg_sur','Dataset-2_with_Mg_sur','Dataset-2_without_Mg_sur'},'RowNames',[Data.Input.SZ_names; "Global"]);
disp(t)
catch
    thx=zeros(8,56)   ;
    for n=1:56
        thx(1:size(Data.Output(n).global_h2o_flux_Tg_Ma,1),n)=[Data.Output(n).global_h2o_flux_Tg_Ma];
    end
    thx(end,:) = sum(thx,2);

    disp('Error - Something did not work.')

    t=array2table(thx([1,2,6],:)','VariableNames',{'Dataset-1_w_Mg_sur','Dataset-2_w_Mg_sur','Dataset-2_wo_Mg_sur'},'RowNames',[Data.Input.SZ_names; "Global"]);
    disp(t)
    disp('Error - Something did not work.')
    disp('Make sure the correct path to the output files is set and all subduction zones results are in the output folder')

end
clearvars -except  Data

