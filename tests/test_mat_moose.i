
       
            
        [Materials]

        [elasticity_tensor]
            type = ComputeIsotropicElasticityTensor
            youngs_modulus = 1e12 # 
            #bulk_modulus = 0.4e6
            poissons_ratio = 1E-6 #  
            compute=true
            block = '2 4 6 8 9 10 11 12 15 16 17'
        []
        
                    [densitybasement]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  2 6 9 10 
                        prop_values = 2600 
                    []
                    

                    [densitysediments]
                        type = GenericConstantMaterial
                        prop_names = 'density'
                        block =  4 8 11 12 
                        prop_values = 2300 
                    []
                    
        [stress]
        type = ComputeLinearElasticStress
        block= '2 4 6 8 9 10 11 12 15 16 17'
        []

        [ini_stress]
              type = ComputeEigenstrainFromInitialStress
              eigenstrain_name = ini_stress
        
            #Initial stress scenario, initial stress tensor is calculated as a factor of vertical gravitational stress
                initial_stress = '1 1 1 1 1 1 1 1 1'
                initial_stress_aux = 'initStressXX initStressXY initStressXZ initStressYX initStressYY initStressYZ initStressZX initStressZY initStressZZ'
                block = '2 4 6 8 9 10 11 12 15 16 17'
        [] 


    
        
        