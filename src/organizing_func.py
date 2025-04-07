import numpy as np
import pandas as pd
import sys

def XYZ_stress_tensor(azimuth_deg, z, rho, g=9.81, k_H=0.6, k_h=0.4, print_output=False):
        """
        Generalized stress tensor as a function of depth, considering faulting regimes
        and aligning the tensor with the azimuth of the maximum principal horizontal stress.
        #### IMPORTANT: THIS OUTPUTS A STRESS TENSOR THAT CONFORMS TO AN [ X Y Z ] COORDINATE SPACE #### 
        #### MANY TEXTBOOK EXAMPLES USE A NorthEastDown (NED) system e.g. Zoback 2010
        
        Parameters:
        z (float): Depth (negative meters below the surface).
        rho (float): Density of the overburden rock (kg/m^3).
        azimuth (float): Azimuth of maximum horizontal stress (degrees clockwise from north), North parallel with Y axis.
        g (float): Gravitational acceleration (m/s^2), default is 9.81 m/s^2. # this is not used in this version, was important for backwards compatibility
        k_H (float): Proportionality constant for maximum horizontal stress (relative to vertical stress).
        k_h (float): Proportionality constant for minimum horizontal stress (relative to vertical stress).
        
        Returns:
        rotated_tensor (numpy.ndarray): 3x3 rotated stress tensor matrix based on the faulting regime in [X Y Z} coordinate space.
        """
        if z >= 0:
            raise ValueError("Depth (z) should be a negative value representing meters below the surface")
        if print_output:
            print("Azi, SHMax: ", azimuth_deg)
        azimuth_deg = azimuth_deg+90
        # Calculate vertical stress due to gravity
        sigma_v = rho * (z)#rho * g * (-z)  # Vertical stress increases with depth
        
        # Calculate horizontal stresses
        sigma_H = k_H * sigma_v   # Maximum horizontal stress
        sigma_h = k_h * sigma_v   # Minimum horizontal stress
        if print_output:
            print("Sigma_H - K : ", sigma_H)
            print("Sigma_h - k : ", sigma_h)
            print("Sigma V: ", sigma_v)
        
        azimuth_rad = np.radians((azimuth_deg))
        # Determine the faulting regime
        if sigma_v < sigma_H < sigma_h:
            # Normal faulting regime: σ1 = σv, σ2 = σH, σ3 = σh
            sigma_xx, sigma_yy, sigma_zz = sigma_v, sigma_H, sigma_h
            alpha, beta, gamma = np.radians(azimuth_deg-90), np.radians(90), np.radians(0)
        elif sigma_H < sigma_v < sigma_h:
            # Strike-slip faulting regime: σ1 = σH, σ2 = σv, σ3 = σh
            sigma_xx, sigma_yy, sigma_zz = sigma_H, sigma_v, sigma_h
            alpha, beta, gamma = azimuth_rad, np.radians(0), np.radians(90)
        elif sigma_H < sigma_h < sigma_v:
            # Reverse faulting regime: σ1 = σH, σ2 = σh, σ3 = σv
            sigma_xx, sigma_yy, sigma_zz = sigma_H, sigma_h, sigma_v
            alpha, beta, gamma = azimuth_rad, np.radians(0), np.radians(0)
        else:
            raise ValueError("Invalid stress regime detected. Check input parameters.")
        
        # Assume no shear stresses (off-diagonal terms = 0)
        sigma_xy = sigma_xz = sigma_yz = 0.0
        
        # Build the initial stress tensor in the global coordinate system
        stress_tensor = np.array([
            [sigma_xx, sigma_xy, sigma_xz],
            [sigma_xy, sigma_yy, sigma_yz],
            [sigma_xz, sigma_yz, sigma_zz]
        ])
        if print_output:
            print("Principle Stress Tensor: \n",stress_tensor.round(3))
        # Convert azimuth to radians
        
        
        # Define the rotation matrix for the azimuth
        rotation_matrix = np.array([
            [np.cos(alpha)*np.cos(beta), np.sin(alpha)*np.cos(beta), -np.sin(beta)],
            [(np.cos(alpha)*np.sin(beta)*np.sin(gamma)) - (np.sin(alpha)*np.cos(gamma)),  (np.sin(alpha)*np.sin(beta)*np.sin(gamma)) + (np.cos(alpha)*np.cos(gamma)), (np.cos(beta)*np.sin(gamma))],
            [(np.cos(alpha)*np.sin(beta)*np.cos(gamma)) + (np.sin(alpha)*np.sin(gamma)),  (np.sin(alpha)*np.sin(beta)*np.cos(gamma)) - (np.cos(alpha)*np.sin(gamma)), (np.cos(beta)*np.cos(gamma))]
        ])
        if print_output:
            print("Rotation Matrix: \n", rotation_matrix.round(3))
        # Rotate the stress tensor to align with the azimuth
        rotated_tensor = rotation_matrix.T @ stress_tensor @ rotation_matrix
        if print_output:
            print("Rotated Matrix: \n",rotated_tensor.round(3))
        return rotated_tensor


    #######    String formatting functions   

def combine_columns_to_string(df: pd.DataFrame, index: int, columns: list) -> str:
            """
            Combine specified columns of a DataFrame row into a comma-separated string.
            Strips any apostrophes and brackets within the values, except for the ones closing the final string.
        
            :param df: pandas DataFrame
            :param index: Row index to retrieve values from
            :param columns: List of column names to include in the output string
            :return: Formatted string with combined values
            """
            values = [str(df.at[index, col]).replace("'", "").replace("[", "").replace("]", "").replace(",","") for col in columns]
            return f"'{', '.join(values).replace(",","")}'"
        
def surf_list(prefix="Surface_", start=1, end=2):
            """
            Return a single quote enclosed sequentionally number list of whatever you pass into the 'prefix' argument

            Parameters:
            prefix (str) : text of list prefix
            start (int) : start number of zero padded 4 digit list
            end (int) : end number... like start
            """
            boundary_list = ' '.join([f'{prefix}{i:04d}' for i in range(start, end+1)])
            #print(boundary_list)
            return boundary_list