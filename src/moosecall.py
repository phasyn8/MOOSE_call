'''MOOSE input file builder class:
    Intended use for distributed probabalistic model simulation or,
    if you would just like to create MOOSE inputs from template simulations'''

"""
Parameters:
    path, mesh_path, meshfile, filename, version, SHmax, Shmin, Sv, surface_elev
    
    **kwargs: [TODO, create different output methods]
    
    kwargs:
        
"""
import os # this could also be used to automatically execute the MOOSE simulations
import pandas as pd
import numpy as np
import sys
sys.path.insert(1, './') # for dependencies
import organizing_func as func

class mooseCall(object):


# ++++++++++++++++++++++++++++++ INIT ++++++++++++++++++++++++++++++++++++++

    def __init__(self, filename, file_path, mesh_path, version, parameters_df, surface_elev, **kwargs):
        #Keyword Args assingment
            
            if kwargs.get("density_dict") == None:
                self.density_dict = {}
            else:
                self.density_dict = kwargs.get("density_dict")

            if kwargs.get("stress_tensor") == None:
                self.stress_tensor = np.eye(3)
            else:
                self.stress_tensor = kwargs.get("stress_tensor")

            if kwargs.get("description") == None:
                self.description= "No description"     
            else:
                self.description = kwargs.get("description") # this is a short description
            
            if kwargs.get('simulation_name') == None:
                self.simulation_name = "No Simulation Name"
            else:
                self.simulation_name = kwargs.get('simulation_name')  # this is the project name
    
        #Required arguments
            self.filename = filename
            self.file_path = file_path

            self.write_path = file_path + filename
            self.mesh_path = mesh_path
            self.version = version
            self.surface_elev = surface_elev

            # Pull Surface and Volume parameters from a mesh dataframe.    
            self.parameters = parameters_df

            vol_list = parameters_df.loc[:, parameters_df.columns.str.contains("Volumes")].columns.tolist()
            #meshLOG[vol_list[0]]  # this prints the volume list
            #vol_list = meshLOG.loc[:, meshLOG.columns.str.contains("Volumes")].columns.tolist()
            #meshLOG[vol_list]

            # init vars 
            self.meshBlockIsBuilt_flag = False
            self.start_blockID = 10000 # starts at an arbitrarily high number to avoid conflict with the fragmented meshblocks
            self.current_blockID = self.start_blockID 

            #====================Surface and Fault Blocks ===========
                    # these are abstracted to avoid conflicts with other functions
            meshLOG = self.parameters  
            current_blockID = self.start_blockID
            fault_meshBlocks = [] # This needs to be populated from the meshblock builder
            surface_meshBlocks = [] # This list needs to be populated from the meshblock builder

            fault_list = meshLOG.loc[:, meshLOG.columns.str.contains("Fault_Elements")].columns.tolist()
            fault_names = [item.split(' ')[0] for item in fault_list]

            surface_list = meshLOG.loc[:, meshLOG.columns.str.contains("Surface_Elements")].columns.tolist()
            surface_names = [item.split(' ')[0] for item in surface_list]


            """  SURFACE LISTS """
            #for surface_name, surface_sidesets, in zip(surface_names, meshLOG[surface_list].values.tolist()):
            for idx, surface_name in enumerate(surface_names):
                #surface_sidesets = func.combine_columns_to_string(meshLOG, index=0, columns=[surface_list[idx]]) # this is only useful for mesh blocks
                surface_meshBlocks.append(f"{surface_name}_surface")
                current_blockID = current_blockID+1


            """ FAULT LISTS """
            #for fault_name, fault_sidesets, in zip(fault_names, meshLOG[fault_list].values.tolist()):
            for idx, fault_name in enumerate(fault_names):
                #fault_sidesets = func.combine_columns_to_string(meshLOG, index=0, columns=[fault_list[idx]]) # this is only useful for mesh blocks
                fault_meshBlocks.append(f"{fault_name}_fault")
                current_blockID = current_blockID+1
                #====================END MESHBLOCK ACCOUNTING ===========
            
            # Building strings for MOOSE Block args
            self.vol_meshBlocks = func.combine_columns_to_string(df=parameters_df, index=0, columns=vol_list)
            self.surface_meshBlocks = f"'{' '.join(surface_meshBlocks)}'" # joining list into sting to pass as MOOSE argument 
            self.fault_meshBlocks = f"'{' '.join(fault_meshBlocks)}'"

            meshBlocks_list = [self.vol_meshBlocks[1:-1], self.surface_meshBlocks[1:-1], self.fault_meshBlocks[1:-1]]

            self.all_meshBlocks = f"'{' '.join(meshBlocks_list)}'"

            azimuth_deg = parameters_df['Azimuth SHMax'].values[0]
            k_h = parameters_df['Sv ratio k_Shmin'].values[0]
            k_H = parameters_df['Sv ratio K_SHMax'].values[0]

            self.stress_tensor = func.XYZ_stress_tensor(azimuth_deg=azimuth_deg, z=-1, rho=1, g=9.81, k_H=k_H, k_h=k_h, print_output=False)


                #use this to parse a dictionary input...
            #formation_map = {key: value[0] for key, value in formation_map_type.items()}\

    def obj_parameters(self):
         print("Write path : ", self.write_path)
         print("Mesh path : ", self.mesh_path)
         print("Version : ", self.version)
         print("Surface Elevation : ", self.surface_elev)
         print("Formation Densities dictionary:", self.density_dict)

         print("All Mesh Blocks : ", self.all_meshBlocks)
         print("Surface Mesh Blocks : ", self.surface_meshBlocks)
         print("Fault Mesh Blocks : ", self.fault_meshBlocks)
         
         print("Stress Tensor :", self.stress_tensor)
         
         print("Parameter Dataframe:")
         return self.parameters

    

#+++++++++++++++++++++++ THIS IS THE START OF THE SIMULATION TEMPLATES +++++++++++++++++++++++++++
    
    def preGravity_stress(self):
        """
        {self.headerBlock}
        {self.globalBlock}
        {self.meshBlock}
        {self.varsBlock}
        {self.auxVarsBlock}
        {self.physicsBlock}
        {self.kernelsBlock}
        {self.auxKernelsBlock}  
        {self.bcBlock}
        {self.functionBlock}
        {self.icBlock}
        {self.materialBlock}
        {self.userobjectBlock}
        {self.postBlock}
        {self.preconditionBlock}
        {self.execBlock}
        {self.outputBlock}
        """
        suffix = "_pre_grav"
        # Simulation blocks construction
        self.general_headers_Blk()
        self.general_globals_Blk()
        if self.meshBlockIsBuilt_flag:
            pass
        else:
            self.generalized_mesher_mesh_Blk()
        self.general_varaibles_Blk()
        self.preGravity_AuxVars_Blk()
        self.solidMech_quasiStatic_physics_Blk()
        self.preGravity_Kernels_Blk()
        self.preGravity_AuxKernels_Blk()
        self.atmos_geomech_BC_Blk()
        self.no_functions_Blk()
        self.no_ic_Blk()
        self.preGravity_Projected_Stress_materials_Blk(density_dict=self.density_dict)
        self.no_UO_Blk()
        self.no_PostProcessor()
        self.general_precond_Blk()
        self.transient_exec_Blk()
        self.output_Block(exodus=True, csv=False, csv_deliminator=',')
        
        # combine and write MOOSE input
        self.moosebuilder(suffix=suffix)
        self.preGravity_exodus = f'{self.write_path[:-2]}{suffix}_out.e'
        print(f"Exported pre gravity calc{self.write_path[:-2]}{suffix}.i")


    def fault_stability_projected_stress_sim(self):
        
        """
        {self.headerBlock}
        {self.globalBlock}
        {self.meshBlock}
        {self.varsBlock}
        {self.auxVarsBlock}
        {self.physicsBlock}
        {self.kernelsBlock}
        {self.auxKernelsBlock}  
        {self.bcBlock}
        {self.functionBlock}
        {self.icBlock}
        {self.materialBlock}
        {self.userobjectBlock}
        {self.postBlock}
        {self.preconditionBlock}
        {self.execBlock}
        {self.outputBlock}
        """
        suffix='_fault_stability'
        #general_headers_Blk(self)
        # Simulation blocks construction
        self.general_headers_Blk()
        self.general_globals_Blk()
        if self.meshBlockIsBuilt_flag:
            pass
        else:
            self.generalized_mesher_mesh_Blk()
        self.general_varaibles_Blk()
        self.postGravity_AuxVars_Blk()
        self.solidMech_quasiStatic_physics_Blk()
        self.postGravity_Kernels_Blk()
        self.postGravity_AuxKernels_Blk()
        self.atmos_geomech_BC_Blk()
        self.no_functions_Blk()
        self.ic_initStress_Blk() #self.no_ic_Blk()
        self.postGravity_Projected_stress_materials_Blk()
        self.solution_UO_Blk(solution_exodus=self.preGravity_exodus)#self.no_UO_Blk()
        self.ts_td_sf_PostProcessors()#self.no_PostProcessor()
        self.general_precond_Blk()
        self.steady_exec_Blk()
        self.output_Block(exodus=True, csv=True, csv_deliminator=',')
        
        # combine and write MOOSE input
        self.moosebuilder(suffix=suffix)

        print(f"Exported Post gravity calc projected stress fault stability to : {self.write_path}")





#+++++++++++++++++++++++   THIS IS THE START OF THE BLOCK BUILDING TEMPLATES    +++++++++++++++++++++++++++

    def moosebuilder(self, suffix='no_suffix'):
        
        """
        Combines the input blocks into an input file 
        """
        
        # Input Builder

        MOOSE_input = f"""
        {self.headerBlock}
        {self.globalBlock}
        {self.meshBlock}
        {self.varsBlock}
        {self.auxVarsBlock}
        {self.physicsBlock}
        {self.kernelsBlock}
        {self.auxKernelsBlock}  
        {self.bcBlock}
        {self.functionBlock}
        {self.materialBlock}
        {self.icBlock}
        {self.userobjectBlock}
        {self.postBlock}
        {self.preconditionBlock}
        {self.execBlock}
        {self.outputBlock}
        """
        f = open(f'{self.write_path[:-2]}{suffix}.i', 'w+')
            
        f.write(MOOSE_input)
        f.close()
            
        print(f"Successfully exported MOOSE input file {self.filename} to {self.file_path}")

#__________________________________________________________ HEADERS BLOCK _________________________________________________________________
    def general_headers_Blk(self):
         
        self.headerBlock= f"""
            #{self.simulation_name}
            #{self.description}
            
            #Version is {self.version}
            #Mesh Path : {self.mesh_path}
            
            #Direction SHmax (AZI Clockwise from North) = {self.parameters['Azimuth SHMax'].values}
            #Sv / SHMax ration = {self.parameters['Sv ratio K_SHMax'].values}
            #Sv / Shmin ration = {self.parameters['Sv ratio k_Shmin'].values}

            #Surface elevation = {self.surface_elev}
        
        """
#__________________________________________________________ GLOBALS BLOCK _________________________________________________________________
    def general_globals_Blk(self):
        
        self.globalBlock = f"""
        # Instantiate Displacement vars
        [GlobalParams]
            displacements = 'disp_x disp_y disp_z'
        []
        """
#__________________________________________________________ MESH BLOCK _________________________________________________________________

    def generalized_mesher_mesh_Blk(self):
        
        meshLOG = self.parameters


        fault_list = meshLOG.loc[:, meshLOG.columns.str.contains("Fault_Elements")].columns.tolist()
        fault_names = [item.split(' ')[0] for item in fault_list]

        surface_list = meshLOG.loc[:, meshLOG.columns.str.contains("Surface_Elements")].columns.tolist()
        surface_names = [item.split(' ')[0] for item in surface_list]

        faults_lowerDBlock_list = []
        surface_lowerDBlock_list = []

        current_msh = 'msh' # this instantiates the first mesh name and will be updated for each modified mesh through the function


        """  SURFACE LISTS """
        #for surface_name, surface_sidesets, in zip(surface_names, meshLOG[surface_list].values.tolist()):
        for idx, surface_name in enumerate(surface_names):
            surface_sidesets = func.combine_columns_to_string(meshLOG, index=0, columns=[surface_list[idx]])
            next_msh = f'{surface_name}Surface'
            surface_lowerDBlock_list.append(f"""
                    [{next_msh}]
                        input={current_msh}
                        type = LowerDBlockFromSidesetGenerator
                        new_block_id = {self.current_blockID}
                        new_block_name = "{surface_name}_surface"
                        sidesets = {surface_sidesets}
                    []"""
            )
            #self.surface_meshBlocks.append(f"{surface_name}_surface")
            current_msh = f'{surface_name}Surface'
            self.current_blockID = self.current_blockID+1


        """ FAULT LISTS """
        #for fault_name, fault_sidesets, in zip(fault_names, meshLOG[fault_list].values.tolist()):
        for idx, fault_name in enumerate(fault_names):
            fault_sidesets = func.combine_columns_to_string(meshLOG, index=0, columns=[fault_list[idx]])
            next_msh = f'{fault_name}Fault'
            faults_lowerDBlock_list.append(f"""
                    [{next_msh}]
                        input={current_msh}
                        type = LowerDBlockFromSidesetGenerator
                        new_block_id = {self.current_blockID}
                        new_block_name = "{fault_name}_fault"
                        sidesets = {fault_sidesets}
                    []"""
            )
            #self.fault_meshBlocks.append(f"{fault_name}_fault")
            current_msh = f'{fault_name}Fault'
            self.current_blockID = self.current_blockID+1
            
            
        fault_lower_dimensionalBlocks = "\n".join(faults_lowerDBlock_list) # converts the list to a formatted string
        surface_lower_dimensionalBlocks = "\n".join(surface_lowerDBlock_list) # converts the list to a formatted string                
            
        self.meshBlock = f"""
        [Mesh]
            [msh]
            type = FileMeshGenerator
            file = '{self.mesh_path}'
            #construct_side_list_from_node_list=false # FOR MESHES WITHOUT IMPLICITLY DEFINED SIDESETS
            []
            {surface_lower_dimensionalBlocks}
            {fault_lower_dimensionalBlocks} 
        []   #  ------------ END OF MESH BLOCK
        """
        self.meshBlockIsBuilt_flag = True
    
    
    def cal_create_tracked_mesh_block(self, dim2_elements_lists, mesh_path, meshfile, filename):

        """
        Creating MOOSE mesh blocks with a list of dimention 2 lists from a fragmenting mesh building operation

        This is specific to the Californie
    
        """

        beld_list = str(dim2_elements_lists[0])[1:-1]
        teg_list = str(dim2_elements_lists[1])[1:-1]
        dulk_list = str(dim2_elements_lists[2])[1:-1]
        vier_list = str(dim2_elements_lists[3])[1:-1]
        zee_list = str(dim2_elements_lists[4])[1:-1]
        namu_list = str(dim2_elements_lists[5])[1:-1]
        nsg_list = str(dim2_elements_lists[6])[1:-1]

        meshBlock = f"""
        [Mesh]
            [msh]
            type = FileMeshGenerator
            file = '{mesh_path}{meshfile}'
            #construct_side_list_from_node_list=false # FOR MESHES WITHOUT IMPLICITLY DEFINED SIDESETS
            []

        # MESH GENERATOR NEEDS CASCADING MESHES TO COMBINE MULTIPLE LowerDimentionalBlocks i.e, one solve feeds into the next
            [TegelenFault]
                input=msh
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 5
                new_block_name = "Tegelen_fault"
                sidesets = {teg_list}
            []
            [DulkenerFault]
                input=TegelenFault
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 6
                new_block_name = "Dulkener_fault"
                sidesets= {dulk_list}
            []
            [BelfeldFault]
                input=DulkenerFault
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 7
                new_block_name = "Belfeld_fault"
                sidesets= {beld_list}
            []
            [ViersenFault]
                input=BelfeldFault
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 8
                new_block_name = "Viersen_fault"
                sidesets= {vier_list}
            []
            [ZeelandFM]
                input=ViersenFault
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 9
                new_block_name = "Zeeland_Fm"
                sidesets= {zee_list}
            []
            [NamurianFm]
                input=ZeelandFM
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 10
                new_block_name = "ZeelandFm"
                sidesets= {namu_list}
            []  
            [NorthSeaGroup]
                input=NamurianFm
                type = LowerDBlockFromSidesetGenerator
                new_block_id = 11
                new_block_name = "NorthSeaGroup"
                sidesets= {nsg_list}
            []
            [] # End Mesh Block ==================================================================================================================

        """
        self.meshBlock = f"!include {mesh_path+filename}"
        f = open(mesh_path+filename, 'w+')
        f.write(meshBlock)
        f.close()
        self.meshBlockIsBuilt_flag = True

#__________________________________________________________ VARIABLES BLOCK _________________________________________________________________
    def general_varaibles_Blk(self):
        self.varsBlock = f"""
            [Variables]
                [disp_x]
                []
                [disp_y]
                []
                [disp_z]
                []
            [] # End Variables Block ============================================================================================================
            """
#__________________________________________________________ AuxVARIABLES BLOCK _________________________________________________________________

    def preGravity_AuxVars_Blk(self):
        self.auxVarsBlock = f"""

        [AuxVariables]

            # Single Stress tensor components for gravitational loading only   

            [stress_zz]
                order = CONSTANT
                family = MONOMIAL
            []

                # Initial Stress Vars computed from Stress_ZZ due to gravitational loading

                [initStressX_X]
                    order = CONSTANT
                    family = MONOMIAL
                []
                [initStressX_Y]
                    order = CONSTANT
                    family = MONOMIAL
                []
                [initStressX_Z]
                    order = CONSTANT
                    family = MONOMIAL
                []

                    [initStressY_X]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressY_Y]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressY_Z]
                        order = CONSTANT
                        family = MONOMIAL
                    []

                        [initStressZ_X]
                            order = CONSTANT
                            family = MONOMIAL
                        []
                        [initStressZ_Y]
                            order = CONSTANT
                            family = MONOMIAL
                        []
                        [initStressZ_Z]
                            order = CONSTANT
                            family = MONOMIAL
                        []
        [] # End AuxVariables Block =====================================================================================================
        """

    def postGravity_AuxVars_Blk(self):
        self.auxVarsBlock = f"""
        [AuxVariables]
            # Stress scalers solved for element directions      
            [normal_vec_magnitude]
                order = CONSTANT
                family = MONOMIAL
            []
            [normal_stress_magnitude]
                order = CONSTANT
                family = MONOMIAL
            []
            [eff_normal_stress_magnitude]
                order = CONSTANT
                family = MONOMIAL
            []
                    [./maxprincipal]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
                    [./midprincipal]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
                    [./minprincipal]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
                    [./hydrostaticRankTwo]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
                    [./hydrostatic_FW]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
                    [./true_ShearStress]
                        order = CONSTANT
                        family = MONOMIAL
                    [../]
            
                            # Element Normal vector components from porous flow module --- TODO: try to make these vector vars! ---
                            [normal_x]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [normal_y]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [normal_z]
                                order = CONSTANT
                                family = MONOMIAL
                            []

            # Normalized normal vector components  --- TODO: try to make these vector vars! ---
            [normal_unit_dir_x]
                order = CONSTANT
                family = MONOMIAL
            []
            [normal_unit_dir_y]
                order = CONSTANT
                family = MONOMIAL
            []
            [normal_unit_dir_z]
                order = CONSTANT
                family = MONOMIAL
            []
            # Tau Strike/Dip vars
            [tau_s]
                order = CONSTANT
                family = MONOMIAL
            []
            [tau_d]
                order = CONSTANT
                family = MONOMIAL
            []

                    # Normal stress vector components --- TODO: try to make these vector vars! ---
                    [normal_stress_vec_x]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [normal_stress_vec_y]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [normal_stress_vec_z]
                        order = CONSTANT
                        family = MONOMIAL
                    []

                            #Stress tensor components       
                            [stress_xx]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_yy]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_zz]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_xy]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_xz]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_yx]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_yz]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_zx]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                            [stress_zy]
                                order = CONSTANT
                                family = MONOMIAL
                            []
                

                    # Initial Stress Vars computed from stress_ZZ due to gravitational loading in preGravity simulation
                    [initStressXX]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressXY]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressXZ]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressYX]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressYY]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressYZ]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressZX]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressZY]
                        order = CONSTANT
                        family = MONOMIAL
                    []
                    [initStressZZ]
                        order = CONSTANT
                        family = MONOMIAL
                    []

                [./MaxShear_stress]
                    order = CONSTANT
                    family = MONOMIAL
                []
                [slip_tendency]
                    order = CONSTANT
                    family = MONOMIAL
                []
                [dilation_tendency]
                    order = CONSTANT
                    family = MONOMIAL
                []
                [fracture_suscept]
                order = CONSTANT
                family = MONOMIAL
                [] 

                 
            #Shear vector vars  ---- TODO: try to make these vector vars! -----
            
            [strike_x]
                order = CONSTANT
                family = MONOMIAL
            []
            [strike_y]
                order = CONSTANT
                family = MONOMIAL
            []
            [strike_z]
                order = CONSTANT
                family = MONOMIAL
            []
            [dip_x]
                order = CONSTANT
                family = MONOMIAL
            []
            [dip_y]
                order = CONSTANT
                family = MONOMIAL
            []
            [dip_z]
                order = CONSTANT
                family = MONOMIAL
            []
            [rake_x]
                order = CONSTANT
                family = MONOMIAL
            []
            [rake_y]
                order = CONSTANT
                family = MONOMIAL
            []
            [rake_z]
                order = CONSTANT
                family = MONOMIAL
            []
            [rake_angle]
                order = CONSTANT
                family = MONOMIAL
            [] 
            [maxShear_angle]
                order = CONSTANT
                family = MONOMIAL
            [] 
            [shear_x]
                order = CONSTANT
                family = MONOMIAL
            []
            [shear_y]
                order = CONSTANT
                family = MONOMIAL
            []
            [shear_z]
                order = CONSTANT
                family = MONOMIAL
            []

        [] # End AuxVariables Block =====================================================================================================
        """
#__________________________________________________________ KERNELS BLOCK _________________________________________________________________
    def preGravity_Kernels_Blk(self):

        self.kernelsBlock = f"""
            [Kernels]
                [gravity]
                    type = Gravity
                    variable = disp_z
                    value = -9.81E-6  # Gravity in MPa
                    block = {self.all_meshBlocks}
                []
            [] # End Kernels Block ===============================================================================================================
        """
    
    def postGravity_Kernels_Blk(self):

        self.kernelsBlock = f"""
            [Kernels]
                #[gravity]
                #    type = Gravity
                #    variable = disp_z
                #    value = -9.81E-6  # Gravity in MPa
                #    block = {self.all_meshBlocks}
                #[]
            [] # End Kernels Block ===============================================================================================================
        """



#__________________________________________________________ AuxKERNELS BLOCK _________________________________________________________________

    def preGravity_AuxKernels_Blk(self):

        K_XX, K_XY, K_XZ, K_YX, K_YY, K_YZ, K_ZX, K_ZY, K_ZZ = self.stress_tensor.ravel()


        self.auxKernelsBlock = f"""
            [AuxKernels]
                # Element specific Normals for each component vector
                [./Stress_zz]
                    type = RankTwoAux
                    rank_two_tensor = stress
                    variable = stress_zz
                    index_i = 2
                    index_j = 2
                    execute_on = TIMESTEP_END
                    block = {self.all_meshBlocks}
                [../]

        #Initial stress function for stiff multi layered material, stress_zz is calculated by gravitational force kernel

        [init_Stress_XX]
                type = ParsedAux
                variable = initStressX_X
                coupled_variables = 'stress_zz'
                expression = 'abs(stress_zz)*{K_XX}' # e.g : K_h = 0.8 Sv
                #execute_on = 'initial timestep_begin'
                block = {self.all_meshBlocks}
        []
        [init_Stress_XY]
                type = ParsedAux
                variable = initStressX_Y
                coupled_variables = 'stress_zz'
                expression = 'abs(stress_zz)*{K_XY}'  
                block = {self.all_meshBlocks}
        []
        [init_Stress_XZ]
                type = ParsedAux
                variable = initStressX_Z
                coupled_variables = 'stress_zz'
                expression = 'abs(stress_zz)*{K_XZ}'  
                block = {self.all_meshBlocks}
        []

                    [init_Stress_YX]
                            type = ParsedAux
                            variable = initStressY_X
                            coupled_variables = 'stress_zz'
                            expression = 'abs(stress_zz)*{K_YX}'  
                            block = {self.all_meshBlocks}
                    []
                    [init_Stress_YY]
                            type = ParsedAux
                            variable = initStressY_Y
                            coupled_variables = 'stress_zz'
                            expression = 'abs(stress_zz)*{K_YY}' 
                            block = {self.all_meshBlocks}
                    []
                    [init_Stress_YZ]
                            type = ParsedAux
                            variable = initStressY_Z
                            coupled_variables = 'stress_zz'
                            expression = 'abs(stress_zz)*{K_YZ}' 
                            block = {self.all_meshBlocks}
                    []

                                [init_Stress_ZX]
                                        type = ParsedAux
                                        variable = initStressZ_X
                                        coupled_variables = 'stress_zz'
                                        expression = 'abs(stress_zz)*{K_ZX}'  # This is Zero for no shear model as one principle stress is the Vertical stress
                                        block = {self.all_meshBlocks}
                                []
                                [init_Stress_ZY]
                                        type = ParsedAux
                                        variable = initStressZ_Y
                                        coupled_variables = 'stress_zz'
                                        expression = 'abs(stress_zz)*{K_ZY}' # No shear here as well
                                        block = {self.all_meshBlocks}
                                []
                                [init_stress_ZZ]
                                        type = FunctorAux
                                        functor = 'stress_zz'
                                        variable = "initStressZ_Z"
                                        # this auxkernel must execute before the variable functor is used
                                        # in the FunctorCoordinatesFunctionAux if we want y to be up to date!
                                        #execute_on = INITIAL # this does not have the expected result
                                []
    
    [] # End AuxKernels Block ============================================================================================================
        """

    def postGravity_AuxKernels_Blk(self):

        self.auxKernelsBlock = f"""
        [AuxKernels]
                # Element specific Normals for each component vector
                [ElemNormal_x]
                    type = PorousFlowElementNormal
                    variable = normal_x
                    component = x
                    #3D_default = '0 0 1'
                    block = {self.all_meshBlocks}
                []
                [ElemeNormal_y]
                    type = PorousFlowElementNormal
                    variable = normal_y
                    component = y
                    #3D_default = '0 0 1'
                    block = {self.all_meshBlocks}
                []
                [ElemNormal_z]
                    type = PorousFlowElementNormal
                    variable = normal_z
                    component = z
                    #3D_default = '0 0 1'
                    block = {self.all_meshBlocks}
                []

        # 9 components for Stress tensor, RankTwo
        [Stress_xx]
            type = RankTwoAux
            rank_two_tensor = stress
            variable = stress_xx
            index_i = 0
            index_j = 0
            execute_on = TIMESTEP_END
            block = {self.all_meshBlocks}
        []
        [Stress_xy]
            type = RankTwoAux
            rank_two_tensor = stress
            variable = stress_xy
            index_i = 0
            index_j = 1
            execute_on = TIMESTEP_END
            block = {self.all_meshBlocks}
        []
        [Stress_xz]
            type = RankTwoAux
            rank_two_tensor = stress
            variable = stress_xz
            index_i = 0
            index_j = 2
            execute_on = TIMESTEP_END
            block = {self.all_meshBlocks}
        []

                [Stress_yy]
                    type = RankTwoAux
                    rank_two_tensor = stress
                    variable = stress_yy
                    index_i = 1
                    index_j = 1
                    execute_on = TIMESTEP_END
                    block = {self.all_meshBlocks}
                []
                [Stress_yx]    
                type = RankTwoAux
                    rank_two_tensor = stress
                    variable = stress_yx
                    index_i = 1
                    index_j = 0
                    execute_on = TIMESTEP_END
                    block = {self.all_meshBlocks}
                []
                [Stress_yz]
                    type = RankTwoAux
                    rank_two_tensor = stress
                    variable = stress_yz
                    index_i = 1
                    index_j = 2
                    execute_on = TIMESTEP_END
                    block = {self.all_meshBlocks}
                []
        
                        [Stress_zx]
                            type = RankTwoAux
                            rank_two_tensor = stress
                            variable = stress_zx
                            index_i = 2
                            index_j = 0
                            execute_on = TIMESTEP_END
                            block = {self.all_meshBlocks}
                        []
                        [Stress_zy]
                            type = RankTwoAux
                            rank_two_tensor = stress
                            variable = stress_zy
                            index_i = 2
                            index_j = 1
                            execute_on = TIMESTEP_END
                            block = {self.all_meshBlocks}
                        []
                        [Stress_zz]
                            type = RankTwoAux
                            rank_two_tensor = stress
                            variable = stress_zz
                            index_i = 2
                            index_j = 2
                            execute_on = TIMESTEP_END
                            block = {self.all_meshBlocks}
                        []

   
    #Shear Stress vector calculations --- there is likely a more straightforward way to address this with a ParsedVectorAux kernels

      [compute_dip_x]
        type = ParsedAux
        variable = dip_x
        coupled_variables = 'normal_unit_dir_y strike_z normal_unit_dir_z strike_y'
        expression = 'normal_unit_dir_y * strike_z - normal_unit_dir_z * strike_y'  # Cross product for dip vector x-component
        block = {self.fault_meshBlocks}
      []

      [compute_dip_y]
        type = ParsedAux
        variable = dip_y
        coupled_variables = 'normal_unit_dir_z strike_x normal_unit_dir_x strike_z'
        expression = 'normal_unit_dir_z * strike_x - normal_unit_dir_x * strike_z'  # Cross product for dip vector y-component
        block = {self.fault_meshBlocks}
      []

      [compute_dip_z]
        type = ParsedAux
        variable = dip_z
        coupled_variables = 'normal_unit_dir_x strike_y normal_unit_dir_y strike_x'
        expression = 'normal_unit_dir_x * strike_y - normal_unit_dir_y * strike_x'  # Cross product for dip vector z-component
        block = {self.fault_meshBlocks}
      []

                [compute_rake_x]
                    type = ParsedAux
                    variable = rake_x
                    coupled_variables = 'maxShear_angle strike_x dip_x'
                    expression = 'cos(maxShear_angle) * strike_x + sin(maxShear_angle) * dip_x'  # Rake x-component
                    block = {self.fault_meshBlocks}
                []

                [compute_rake_y]
                    type = ParsedAux
                    variable = rake_y
                    coupled_variables = 'maxShear_angle strike_y dip_y'
                    expression = 'cos(maxShear_angle) * strike_y + sin(maxShear_angle) * dip_y'  # Rake y-component
                    block = {self.fault_meshBlocks}
                []

                [compute_rake_z]
                    type = ParsedAux
                    variable = rake_z
                    coupled_variables = 'maxShear_angle strike_z dip_z'
                    expression = 'cos(maxShear_angle) * strike_z + sin(maxShear_angle) * dip_z'  # Rake z-component
                    block = {self.fault_meshBlocks}
                []

    [compute_rake_angle]
        type = ParsedAux
        variable = rake_angle
        coupled_variables = 'maxShear_angle'
        #expression = 'atan((maxprincipal - minprincipal)/(minprincipal - maxprincipal))'
        expression = 'maxShear_angle * (180.0/3.141592653589793238463)' # radians to degrees
        block = {self.fault_meshBlocks}
    []
      
                [compute_shear_x]
                    type = ParsedAux
                    variable = shear_x
                    coupled_variables = 'maxShear_angle strike_x'
                    expression = '(cos(maxShear_angle) * strike_x)'
                    block = {self.fault_meshBlocks}
                []

                [compute_shear_y]
                    type = ParsedAux
                    variable = shear_y
                    coupled_variables = 'maxShear_angle strike_y'
                    expression = '(cos(maxShear_angle) * strike_y)'
                    block = {self.fault_meshBlocks}
                []

                [compute_shear_z]
                    type = ParsedAux
                    variable = shear_z
                    coupled_variables = 'maxShear_angle'
                    expression = '(sin(maxShear_angle))'
                    block = {self.fault_meshBlocks}
                []          

    [max_shear_stress_direction]
        type = ParsedAux
        variable = maxShear_angle
        execute_on = TIMESTEP_END
        coupled_variables = 'tau_d tau_s'
        expression = 'atan(tau_d/tau_s)'
        block = {self.fault_meshBlocks}  
    []

            [compute_cross_product_x]
                type = ParsedAux
                variable = strike_x  # The variable to store the cross product result
                coupled_variables = 'normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '(-normal_unit_dir_y)' # Cross product with vertical vector (0, 0, 1)
                block = {self.fault_meshBlocks} 
            []
            [compute_cross_product_y]
                type = ParsedAux
                variable = strike_y  # The variable to store the cross product result
                coupled_variables = 'normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '(normal_unit_dir_x)' # Cross product with vertical vector (0, 0, 1)
                block = {self.fault_meshBlocks} 
            []
            [compute_cross_product_z]
                type = ParsedAux
                variable = strike_z  # The variable to store the cross product result
                coupled_variables = 'normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '0.0' # Cross product with vertical vector (0, 0, 1)
                block = {self.fault_meshBlocks} 
            []

    # ---------------- end max shear vector calc  

                # This should always compute to 1 ---- TODO : DEPRECIATE THIS...---- 
                        [NormalDirectionMagnitude]
                            type = ParsedAux
                            variable = normal_vec_magnitude
                            expression = 'sqrt(normal_x^2 + normal_y^2 + normal_z^2)'
                            coupled_variables = 'normal_x normal_y normal_z'
                            block = {self.fault_meshBlocks} 
                        []
        
    # Aux Variable derivative of the Element normals vector (So they can be used in coupled calculations as the Element specific Normals cannot be directly)
            [NormalDirection_x]
                type = ParsedAux
                variable = normal_unit_dir_x
                expression = 'normal_x / normal_vec_magnitude'
                coupled_variables = 'normal_x normal_vec_magnitude'
                block = {self.fault_meshBlocks} 
            []
        
            [NormalDirection_y]
                type = ParsedAux
                variable = normal_unit_dir_y
                expression = 'normal_y / normal_vec_magnitude'
                coupled_variables = 'normal_y normal_vec_magnitude'
                block = {self.fault_meshBlocks} 
            []
        
            [NormalDirection_z]
                type = ParsedAux
                variable = normal_unit_dir_z
                expression = 'normal_z / normal_vec_magnitude'
                coupled_variables = 'normal_z normal_vec_magnitude'
                block = {self.fault_meshBlocks} 
            []

        [NormalStressVec_x]
                type = ParsedAux
                variable = normal_stress_vec_x
                coupled_variables = 'stress_xx stress_xy stress_xz normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '(normal_unit_dir_x*stress_xx) + (normal_unit_dir_y*stress_xy) + (normal_unit_dir_z*stress_xz)'
                block = {self.fault_meshBlocks} 
        []
        [NormalStressVec_y]
                type = ParsedAux
                variable = normal_stress_vec_y
                coupled_variables = 'stress_yx stress_yy stress_yz normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '(normal_unit_dir_x*stress_yx) + (normal_unit_dir_y*stress_yy) + (normal_unit_dir_z*stress_yz)'
                block = {self.fault_meshBlocks} 
        []
        [NormalStressVec_z]
                type = ParsedAux
                variable = normal_stress_vec_z
                coupled_variables = 'stress_zx stress_zy stress_zz normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z'
                expression = '(normal_unit_dir_x*stress_zx) + (normal_unit_dir_y*stress_zy) + (normal_unit_dir_z*stress_zz)'
                block = {self.fault_meshBlocks} 
        []
        
                    [./TauD] #Scalar dot product of the inverted dip vector (so it points downward) and the normal stress unit vector (NormalForceVec x,y,z)
                            type = ParsedAux
                            variable = tau_d
                            coupled_variables = 'normal_stress_vec_x normal_stress_vec_y normal_stress_vec_z dip_x dip_y dip_z'
                            expression ='(normal_stress_vec_x*(dip_x))+(normal_stress_vec_y*(dip_y))+(normal_stress_vec_z*(dip_z))'
                            block = {self.fault_meshBlocks} 
                    [../]

                    [./TauS] #Scalar dot product of the strike vector and the normal stress unit vector (NormalForceVec x,y,z)
                            type = ParsedAux
                            variable = tau_s
                            coupled_variables = 'normal_stress_vec_x normal_stress_vec_y normal_stress_vec_z strike_x strike_y strike_z'
                            expression ='(normal_stress_vec_x*(strike_x))+(normal_stress_vec_y*(strike_y))+(normal_stress_vec_z*(strike_z))'
                            block = {self.fault_meshBlocks} 
                    [../]

        [./NormalStressMagnitude]
            type = ParsedAux
            variable = normal_stress_magnitude
            block = {self.fault_meshBlocks}  # Only compute this on 2D surfaces (Lower Dimensional Blocks)
            # Calculate the stress magnitude along the normal direction
            expression = 'normal_unit_dir_x*(normal_unit_dir_x*stress_xx + normal_unit_dir_y*stress_xy + normal_unit_dir_z*stress_xz) +
                          normal_unit_dir_y*(normal_unit_dir_x*stress_yx + normal_unit_dir_y*stress_yy + normal_unit_dir_z*stress_yz) +
                          normal_unit_dir_z*(normal_unit_dir_x*stress_zx + normal_unit_dir_y*stress_zy + normal_unit_dir_z*stress_zz)'
            coupled_variables = 'normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz stress_zx stress_yx stress_zy'
        [../]

        [./effectiveNormalStress]
            type = ParsedAux
            variable = eff_normal_stress_magnitude
            block = {self.fault_meshBlocks}  # Only compute this on 2D surfaces (Lower Dimensional Blocks)
            # Calculate the stress magnitude along the normal direction
            expression = 'normal_unit_dir_x*(normal_unit_dir_x*(stress_xx) + normal_unit_dir_y*(stress_xy) + normal_unit_dir_z*(stress_xz)) +
                          normal_unit_dir_y*(normal_unit_dir_x*(stress_yx) + normal_unit_dir_y*(stress_yy) + normal_unit_dir_z*(stress_yz)) +
                          normal_unit_dir_z*(normal_unit_dir_x*(stress_zx) + normal_unit_dir_y*(stress_zy) + normal_unit_dir_z*(stress_zz))+hydrostatic_FW'
            coupled_variables = 'hydrostatic_FW normal_unit_dir_x normal_unit_dir_y normal_unit_dir_z stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz stress_zx stress_yx stress_zy'
        [../]


    # Principle stress calculations plus hydrostatic
        [./maxprincipal]
          type = RankTwoScalarAux
          rank_two_tensor = stress
          variable = maxprincipal
          scalar_type = MaxPrincipal
        [../]
        [./midprincipal]
          type = RankTwoScalarAux
          rank_two_tensor = stress
          variable = midprincipal
          scalar_type = MidPrincipal
        [../]
        [./minprincipal]
          type = RankTwoScalarAux
          rank_two_tensor = stress
          variable = minprincipal
          scalar_type = MinPrincipal
        [../]

        #[./hydrostaticRankTwo] # This does not conform to the hydrostatic pressure gradient, rather hydrostatic developing due to stress
        #  type = RankTwoScalarAux
        #  rank_two_tensor = stress
        #  variable = hydrostaticRankTwo
        #  scalar_type = Hydrostatic
        #[../]

        [hydrostatic_pp_freshwater]
            type = ParsedAux
            variable = hydrostatic_FW
            expression  = '1e-2 * (0-z)' # Approximate hydrostatic in MPa...  must change '0' datum to the model max model elevation
            #coupled_variables = = 'rhoWater' # pressure per meter freshwater
            use_xyzt = true
            block = {self.all_meshBlocks}
        []

        [True_ShearStress]
            type = ParsedAux
            variable = true_ShearStress
            expression = '0.5*((minprincipal+hydrostatic_FW) + (maxprincipal+hydrostatic_FW))'
            coupled_variables = 'minprincipal maxprincipal hydrostatic_FW'
            block = {self.fault_meshBlocks}
        []
            
        # Computes maximum Shear stress on surfaces, though does not output direction of this vector
        [./MaxShearStress] # Need to investigate if this is sensible for element shear stress, would be nice to include direction vector with it
            type = RankTwoScalarAux
            rank_two_tensor = stress
            variable = MaxShear_stress
            scalar_type = MaxShear
            check_boundary_restricted = false # this is true for non-'LowerDimensionalBlocks'
            block = {self.fault_meshBlocks}
        []    

    # Geomechanical calculations on fault surfaces
      [SlipTendency]
        type = ParsedAux
        variable = slip_tendency
        coupled_variables = 'true_ShearStress normal_stress_magnitude hydrostatic_FW'
        expression = 'abs(true_ShearStress)/(abs(normal_stress_magnitude)-abs(hydrostatic_FW))'  # must invert the Maxshear stress to a compressional stress (negative in Moose) beacuase MaxShear outputs absolute value
        #execute_on = timestep_end
        check_boundary_restricted = false
        use_xyzt = true
        block = {self.fault_meshBlocks} 
      [] 

    [DilationTendency]
        type = ParsedAux
        variable = dilation_tendency
        expression = '(minprincipal-normal_stress_magnitude)/(minprincipal-maxprincipal)' # minprinciple = Sigma1 ; maxprinciple = Sigma3
        coupled_variables = 'maxprincipal minprincipal normal_stress_magnitude'
        #execute_on = timestep_end
        check_boundary_restricted = false
        use_xyzt = true
        block = {self.fault_meshBlocks}   
    []

    [Fracture_Suscpetibility]
        type = ParsedAux
        variable = fracture_suscept
        expression = 'abs(eff_normal_stress_magnitude) - abs(true_ShearStress / 0.6)'
        coupled_variables = 'true_ShearStress eff_normal_stress_magnitude'
        check_boundary_restricted = false
        use_xyzt = true
        block = {self.fault_meshBlocks}
    []

    [] # End AuxKernels Block ============================================================================================================
    """

#__________________________________________________________ PHYSICS BLOCK _________________________________________________________________
    def solidMech_quasiStatic_physics_Blk(self):
        self.physicsBlock = f"""
            [Physics]
                [SolidMechanics]
                    [QuasiStatic]
                        [All]
                            strain = SMALL
                            #add_variables = true
                            block = {self.all_meshBlocks}
                            eigenstrain_names = ini_stress
                            displacements = 'disp_x disp_y disp_z'
                            #coupled_variables = 'disp_x, disp_y, disp_z'
                            #generate_output = 'vonmises_stress'
                        []
                    []    
                []
            []  # End Physics Block =============================================================================================================
        """

#__________________________________________________________ BOUNDARY CONDITIONS BLOCK _________________________________________________________________

    def atmos_geomech_BC_Blk(self):

        """
        Creating MOOSE Boundary conditions block with input from a meshLOG entry of mesh boundary elements

        """

        xBC = func.combine_columns_to_string(df=self.parameters, index=0, columns=["West Boundary_Elements","East Boundary_Elements","South Boundary_Elements", "North Boundary_Elements", "Base Boundary_Elements"])
        yBC = func.combine_columns_to_string(df=self.parameters, index=0, columns=["West Boundary_Elements","East Boundary_Elements","South Boundary_Elements", "North Boundary_Elements", "Base Boundary_Elements"])
        zBC = func.combine_columns_to_string(df=self.parameters, index=0, columns=["Base Boundary_Elements"])
        TopBC = func.combine_columns_to_string(df=self.parameters, index=0, columns=["Top Boundary_Elements"])
        
        #xBC = df_picks_index[["West Boundary_Elements","East Boundary_Elements","South Boundary_Elements", "North Boundary_Elements", "Base Boundary_Elements"]]
        #yBC = df_picks_index[["West Boundary_Elements","East Boundary_Elements","South Boundary_Elements", "North Boundary_Elements", "Base Boundary_Elements"]]
        #zBC = df_picks_index[["Base Boundary_Elements"]]
        self.bcBlock = f"""

        #  Must include the boundaries ( N/S , E/W , Top ) - for each X/Y displacement boundary condition, and only Bottom for Z BCs, this when gravity is used for body force
        [BCs]

            [no_x]
            type = DirichletBC
            variable = disp_x
            boundary = {xBC}
            value = 0.0
            []

            [no_y]
            type = DirichletBC
            variable = disp_y
            boundary = {yBC}
            value = 0.0
            []

            [no_z]
            type = DirichletBC
            variable = disp_z
            boundary = {zBC}
            value = 0.0
            []

            [Top_Atmospheric_Pressure]
            type = NeumannBC
            variable = disp_z
            boundary = {TopBC}
            value = -0.101325  # Atmospheric pressure sea level in MPa
            []
        []  # End Boundary conditions Block ==================================================================================================

        """
        #f = open(path+filename, 'w+')
            
        #f.write(bc_Block)
        #f.close()
        #self.bcBlock = f"!include {filename}"

#__________________________________________________________ FUNCTIONS BLOCK ____________________________________________________________
   
    def no_functions_Blk(self):
            
        self.functionBlock = f"""
        [Functions]
        [] # ---------------------- END FUNCTIONS BLOCK -----------------------------
        """
#__________________________________________________________ INITIAL CONDITIONS BLOCK ____________________________________________________________

    def ic_initStress_Blk(self):
        self.icBlock = f"""
        [ICs]
        [ic_initStressX_X]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressXX
            from_variable = 'initStressX_X'
          []
        [ic_initStressX_Y]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressXY
            from_variable = 'initStressX_Y'
        []
        [ic_initStressX_Z]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressXZ
            from_variable = 'initStressX_Z'
          []
        [ic_initStressY_X]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
           variable = initStressYX
            from_variable = 'initStressY_X'
          []
        [ic_initStressY_Y]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressYY
           from_variable = 'initStressY_Y'
          []
        [ic_initStressY_Z]
          type = SolutionIC
         solution_uo = ini_stress_from_gravity_soln
          variable = initStressYZ
          from_variable = 'initStressY_Z'
        []
        [ic_initStressZ_X]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressZX
            from_variable = 'initStressZ_X'
         []
          [ic_initStressZ_Y]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressZY
            from_variable = 'initStressZ_Y'
          []
          [ic_initStressZ_Z]
            type = SolutionIC
            solution_uo = ini_stress_from_gravity_soln
            variable = initStressZZ
            from_variable = 'stress_zz'
          []
        []  # End initial conditions Block ==================================================================================================
        """

    def no_ic_Blk(self):
         
         self.icBlock = f"""
        [ICs]
        []  # End initial conditions Block ==================================================================================================
        """
         
#__________________________________________________________ MATERIALS BLOCK ____________________________________________________________
        
    def postGravity_Projected_stress_materials_Blk(self, density_dict={}):
        
        density_list =[]
        for key, value, in density_dict.items():
            density_list.append(f"""
                    [density{key}]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {value[1]} 
                        prop_values = {value[0]} 
                    []
                    """)
        density_material_properties = "\n".join(density_list)
        
        surface_material_properties = f"""
                    [densitySurface]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {self.surface_meshBlocks} 
                        prop_values = 1
                    []
                    """
        fault_material_properties = f"""
                    [densityfault]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {self.fault_meshBlocks} 
                        prop_values = 1
                    []
                    """
        
        
        self.materialBlock = f"""    
        [Materials]

        [elasticity_tensor]
            type = ComputeIsotropicElasticityTensor
            youngs_modulus = 1e12 # 
            #bulk_modulus = 0.4e6
            poissons_ratio = 1E-6 #  
            compute=true
            block = {self.all_meshBlocks} 
        []
        {density_material_properties}
        {surface_material_properties}
        {fault_material_properties} 
        [stress]
            type = ComputeLinearElasticStress
            block= {self.all_meshBlocks} 
        []

        [ini_stress]
            type = ComputeEigenstrainFromInitialStress
            eigenstrain_name = ini_stress
        
            #Initial stress equilibrium, initial stress tensor is calculated as a factor of vertical gravitational stress and projected onto a non displaced mesh
                initial_stress = '1 1 1 1 1 1 1 1 1'
                initial_stress_aux = 'initStressXX initStressXY initStressXZ initStressYX initStressYY initStressYZ initStressZX initStressZY initStressZZ'
                block = {self.all_meshBlocks}
        [] 
    [] # ______________END MATERIALS BLOCK
        """
    
    
    def preGravity_Projected_Stress_materials_Blk(self, density_dict={}):
        density_list =[]
        for key, value, in density_dict.items():
            density_list.append(f"""
                    [density{key}]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {value[1]} 
                        prop_values = {value[0]} 
                    []
                    """)
        density_material_properties = "\n".join(density_list)
        
        surface_material_properties = f"""
                    [densitySurface]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {self.surface_meshBlocks} 
                        prop_values = 1
                    []
                    """
        fault_material_properties = f"""
                    [densityfault]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  {self.fault_meshBlocks} 
                        prop_values = 1
                    []
                    """

        
        self.materialBlock = f"""    
        [Materials]

        [elasticity_tensor]
            type = ComputeIsotropicElasticityTensor
            youngs_modulus = 1e12 # 
            #bulk_modulus = 0.4e6
            poissons_ratio = 1E-6 #  
            compute=true
            block = {self.all_meshBlocks}
        []

        {density_material_properties}
        {surface_material_properties}
        {fault_material_properties} 

        [stress]
        type = ComputeLinearElasticStress
        block= {self.all_meshBlocks} 
        []

        [ini_stress]
              type = ComputeEigenstrainFromInitialStress
              eigenstrain_name = ini_stress
        
            #Initial stress scenario: stress tensor is calculated as a factor of vertical gravitational stress, will then be the initial condition of post gravity materials block
                initial_stress = '0 0 0 0 0 0 0 0 0'
                #initial_stress_aux = 'initStressXX initStressXY initStressXZ initStressYX initStressYY initStressYZ initStressZX initStressZY initStressZZ'
                block = {self.all_meshBlocks} 
        [] 

    [] #___________________END Materials Block
    """

#__________________________________________________________ USEROBJECTS BLOCK ____________________________________________________________


    def solution_UO_Blk(self, solution_exodus=''):
        
        self.userobjectBlock =f"""
        [UserObjects]
        [ini_stress_from_gravity_soln]
          type = SolutionUserObject
          mesh = '{solution_exodus}'
          system_variables = 'initStressX_X initStressX_Y initStressX_Z initStressY_X initStressY_Y initStressY_Z initStressZ_X initStressZ_Y stress_zz'
          timestep = LATEST
        []
        [] #____________________ END USEROBJECTS
        """

    def no_UO_Blk(self):

        self.userobjectBlock=f"""
        [UserObjects]
        [] #____________________ END UserObjects
        """
#__________________________________________________________ POSTPROCESSOR BLOCK ____________________________________________________________

    def ts_td_sf_PostProcessors(self):
         
        self.postBlock = f"""

        [Postprocessors]
      
        [] #_____________________ END PostProcessors


        [VectorPostprocessors]
        [Ts_value_sampler]
          type = ElementValueSampler
          variable = 'slip_tendency'
          sort_by = id
          block = {self.fault_meshBlocks}
          execute_on = 'TIMESTEP_END'
        []
        [Td_value_sampler]
         type = ElementValueSampler
         variable = 'dilation_tendency'
         sort_by = id
         block = {self.fault_meshBlocks}
         execute_on = 'TIMESTEP_END'
        []
        [Sf_value_sampler]
          type = ElementValueSampler
          variable = 'fracture_suscept'
          sort_by = id
          block = {self.fault_meshBlocks}
          execute_on = 'TIMESTEP_END'
        []
        

        [] #_____________________END VectorPostProcessor

        """

    def no_PostProcessor(self):
         
        self.postBlock = f"""
        [Postprocessors]
      
        [] #_____________________ END PostProcessors
        """


#__________________________________________________________ PRECONDITIONING BLOCK ____________________________________________________________

    def general_precond_Blk(self):
         
        self.preconditionBlock = f"""
        [Preconditioning]
            [SMP]
            type = SMP
            full = true
            []
        [] #____________________ END PRECONDITIONING
        """

#__________________________________________________________ EXECUTION BLOCK _________________________________________________________________
    def steady_exec_Blk(self):

            self.execBlock = f"""
            [Executioner]
            type = Steady
        
            solve_type = 'NEWTON'

            # overshoot control:
            petsc_options_iname = '-pc_type -snes_linesearch_damping'
            petsc_options_value = ' lu .001'
        
            line_search = cp
        
            nl_abs_tol = 1e-6
            nl_rel_tol = 1e-4
        
            l_max_its = 200
            nl_max_its = 200
        
            [] #____________________ END EXECUTION
            """


    def transient_exec_Blk(self):
            self.execBlock = f"""
            [Executioner]
            #type = Steady
            type = Transient
            solve_type = 'NEWTON'

            # overshoot control:
            petsc_options_iname = '-pc_type -snes_linesearch_damping'
            petsc_options_value = ' lu .001'
        
            line_search = cp
        
            nl_abs_tol = 1e-6
            nl_rel_tol = 1e-4
        
            l_max_its = 200
            nl_max_its = 200
            num_steps = 2
            [] #____________________ END EXECUTION
            """



#__________________________________________________________ OUTPUT BLOCK _________________________________________________________________
    def output_Block(self, exodus=True, csv=False, csv_deliminator=','):
          

        _csv = f"""
                [csv]
                type = CSV
                #output_user_objects = true  # Enable exporting UserObjects
                file_base = {self.version} # The base name for the output file
                delimiter = {csv_deliminator}  # Set delimiter to comma for CSV
                []
                """
        _exodus= f"""
                [out]
                type = Exodus
                #elemental_as_nodal = true
                #exodus = true
                []
                """
        
        if exodus & csv:
            out_Block = f"""
            [Outputs]
            {_exodus}
            {_csv}
            []
            """
        elif exodus:
            out_Block = f"""
            [Outputs]
            {_exodus}
            []
            """
        elif csv:
            out_Block = f"""
            [Outputs]
            {_csv}
            []
            """
        self.outputBlock = out_Block