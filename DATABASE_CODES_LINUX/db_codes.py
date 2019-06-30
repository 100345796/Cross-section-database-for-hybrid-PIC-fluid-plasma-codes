def get_collision_info(root, Type, specie1, specie2, process, file_name):
    
    """
    This function reads the information stored in the database of a collision reaction specified by the user.
    That is, it reads a specific file inside the entire databse specified by the user
            
    INPUTS REQUIRED:
        
        -Database_root: the root directory where the database is stored (in the
                        desktop, in a specific folder, in a USB...)
        -Type: type of the desired interaction. With the current state of the database
               this field must be one of: 'BINARY' or 'UNARY'
        -specie1: first of the species involved in the interaction
        -specie2: second specie involved in the interaction
                        
        NOTE_1: the species are introduced with their sign an their electronic state.
                Some examples are: 'e-' for electron
                                   'Xe' for neutral Xe
                                   'Xe+2' for double ionized Xenon--> Symbol+Charge+#charges
                                   'Xe[5s2]' for excited Xe --> Symbol[electronic state of the las valence e-]
                                                   
        NOTE_2: the specied don't necessarily need to be entered in order.
                        
        NOTE_3: If the type UNARY is introduced, the secondary specie should be left as ''
                                                         
        -process: desired reaction. 'i1' for one ionization of 1 e-
                                    'elastic' for elastic collision
                                    'e[...]' for excitation to the [...] electronic configuration
                                    
        -file_name: name of the file that wants to be read (it must be a .txt)
        
        NOTE_3: If the name of the .txt is unknown, leave it as '' and the code
                will automatically read the first file associated with a database or,
                in the case where only one file is stored, it will read this file.
    
    OUTPUT:
        
        -The dictionary structure containing all the information inside the specified file
    """
    
    import yaml
    
    # It must be checked that the order of the variables introduced in the search of the file follow
    # the rules stablish for the database storae system. 
    
    if (specie1 == 'e-' or specie2 == 'e-'):
        if specie1 == 'e-':
            Specie1 = specie1
            Specie2 = specie2
        elif specie2 == 'e-':
            Specie1 = specie2
            Specie2 = specie1
    elif (specie1 == '' or specie2 == ''):
        if specie1 == '':
            Specie1 = specie2
            Specie2 = specie2
        elif specie2 == '':
            Specie1 = specie1
            Specie2 = specie2
    else:
        species = [specie1, specie2]
        species.sort()
        Specie1 = species[0]
        Specie2 = species[1]
        
    # The file root must be adapted for bnary collisions or unitary reactions
   
    if (Type == 'BINARY'):
        if file_name == '':
            file_root = _select_fileroot_filename_unknown(root, Type, Specie1, Specie2, process)
        else:
            file_root = root+"/"+Type+"/"+Specie1+"/"+Specie2+"/"+process+"/"+file_name
    elif (Type == 'UNARY'):
        if file_name == '':
            from _select_fileroot_filename_unknown import select_fileroot_filename_unknown
            file_root = select_fileroot_filename_unknown(root, Type, Specie1, Specie2, process)
        else:
            file_root = root+"/"+Type+"/"+Specie1+"/"+process+"/"+file_name
        
    # The file is now read and stored as dictionary variable:
    
    file = open(file_root,'r')
    info = yaml.load(file)
    
    return(info)
##########################################################################################################
##########################################################################################################
def interpolate_at_velocity(Database_root, Type, specie1, specie2, process, file_name, Velocity):
    """
    This interpolates to obtain the cross section associated with a specific reaction
    at a specific energy.
            
    INPUTS REQUIRED:
        
        -Database_root: the root directory where the database is stored (in the
                        desktop, in a specific folder, in a USB...)
        -Type: type of the desired interaction. With the current state of the database
               this field must be one of: 'BINARY' or 'UNARY'
        -specie1: first of the species involved in the interaction
        -specie2: second specie involved in the interaction
                        
        NOTE_1: the species are introduced with their sign an their electronic state.
                Some examples are: 'e-' for electron
                                   'Xe' for neutral Xe
                                   'Xe+2' for double ionized Xenon--> Symbol+Charge+#charges
                                   'Xe[5s2]' for excited Xe --> Symbol[electronic state of the las valence e-]
                                                   
        NOTE_2: the specied don't necessarily need to be entered in order.
                        
        NOTE_3: If the type UNARY is introduced, the secondary specie should be left as ''
                                                         
        -process: desired reaction. 'i1' for one ionization of 1 e-
                                    'elastic' for elastic collision
                                    'e[...]' for excitation to the [...] electronic configuration
                                    
        -file_name: name of the file that wants to be read (it must be a .txt)
        
        NOTE_3: If the name of the .txt is unknown, leave it as '' and the code
                will automatically read the first file associated with a database or,
                in the case where only one file is stored, it will read this file.
                
        -Velocity: relative velocity between the particles involved. Provided in [m/s]
    
    OUTPUT:
        
        -The interpolate cross-section in [m^2]
    """
    # Dictionary with collision data
    
    collision_info = get_collision_info(Database_root, Type, specie1, specie2, process, file_name)
    
    # Extraction of energy and differential cross section values from dictionary
    
    [Energy_array, CS_array] = _get_E_and_CS_arrays_from_info(collision_info) #[eV, m^2]
    
    # The velocity must be translated to energy to interpolate. The mass of the
    # electron mass unless electrons are not present.
    if specie1 == 'e-' or specie2 == 'e-':
        Specie1 = 'e-'    
    else:
        Specie1 = specie1[0:2]
    
    file_root = Database_root+"/"+"UNARY"+"/"+Specie1+"/"+"specie_info.txt"
    file = open(file_root,'r')
    import yaml
    info = yaml.load(file)
    
    m = info['MASS'] 
        
    E_to_interpolate =0.5*m*Velocity**2 # [J]
    K = 1.6e-19 # [J/eV] Boltzmann constant
    E_to_interpolate = E_to_interpolate/K # [eV]
    
    import numpy as np
    
    # As the cross section data is only provided inside a small range, it must be interpolated outside it:
    if E_to_interpolate<Energy_array.min():
        if specie1 == 'e-' or specie2 == 'e-':
            if process=='elastic':
                f_cs = CS_array[1] + (CS_array[0]-CS_array[1])*(E_to_interpolate - Energy_array[1])/(Energy_array[0]-Energy_array[1]) # elastic dont have threshold
                if f_cs<=0:
                    interpolated_cross_section=0 # it cannot be negtive
                else:
                    interpolated_cross_section=f_cs 
            else:
                interpolated_cross_section = 0 # below the threshold energy collisions dont take place
        else:
            f_cs = CS_array[1] + (CS_array[0]-CS_array[1])*(E_to_interpolate - Energy_array[1])/(Energy_array[0]-Energy_array[1]) # ch_e dont have threshold
            if f_cs<=0:
                interpolated_cross_section=0 # it cannot be negtive
            else:
                interpolated_cross_section=f_cs   
    elif E_to_interpolate>=Energy_array.min() and E_to_interpolate<=Energy_array.max():
        interpolated_cross_section = np.interp(E_to_interpolate, Energy_array, CS_array) # Inside the region, the cross section is interpolated
    elif E_to_interpolate>Energy_array.max():
        f_cs = CS_array[-2] + (CS_array[-1]-CS_array[-2])*(E_to_interpolate - Energy_array[-2])/(Energy_array[-1]-Energy_array[-2]) # outside the range, the cross section is linearly extrapolated
        if f_cs<=0:
            interpolated_cross_section=0 # it cannot be negtive
        else:
            interpolated_cross_section=f_cs   
    
    return interpolated_cross_section

#############################################################################################################
#############################################################################################################
    
def collision_rate_Maxwellian(Database_root, Type, specie1, specie2, process, file_name, Velocity):
    
    """
    This codes computes de collision rate assuming an isotropic system with  
    a given electron flow velocity.
            
    INPUTS REQUIRED:
    
    -Database_root: the root directory where the database is stored (in the
                    desktop, in a specific folder, in a USB...)
    -Type: type of the desired interaction. With the current state of the database
          this field must be one of: 'BINARY' or 'UNARY'
    -specie1: first of the species involved in the interaction
    -specie2: second specie involved in the interaction
        
        NOTE_1: the species are introduced with their sign an their electronic state.
                Some examples are: 'e-' for electron
                                   'Xe' for neutral Xe
                                   'Xe+2' for doble ionized Xenon--> Symbol+Charge+#charges
                                   'Xe[5s2]' for excited Xe --> Symbol[electronic state of the las valence e-]
                                   
        NOTE_2: the specied dont necessarily need to be entered in order.
        
        NOTE_3: If the type UNARY is introduced, the secondary specie should be left as ''
                                         
    -process: desired reaction. 'i1' for one ionization of 1 e-
                                'elastic' for elastic collision
                                'e[...]' for excitation to the [...] electronic configuration
                                
    -file_name: name of the file that wants to be read (it must be a .txt)
    
        NOTE_3: If the name of the .txt is unknown, leave it as '' and the code
                will automatically read the first file associated with a database or,
                in the case where only one file is stored, it will take read this file.
                
    -Velocity: velocity of the flow of electrons. In [eV]
    
    OUTPUT:
    
    -The collision rate at the specified energy. In [m^3/s]
    """
    
    # Dictionary with collision data
    
    collision_info = get_collision_info(Database_root, Type, specie1, specie2, process, file_name)
    
    # Extraction of energy and differential cross section values from dictionary

    [Energy_array, CS_array] = _get_E_and_CS_arrays_from_info(collision_info) #[eV, m^2]
    
    # Important parameteres for the Maxwellian distribution:
    
    m = 9.10938356e-31              # [Kg] electron mass
    K = 1.6e-19                     # [J/eV] Boltzmann constant
    
    # Adjusting of the flow velocity to S.I.:
    
    T = K*Velocity                  # [J] Flow temperature
    Energy_array_J = K*Energy_array            # [J] Energies associated with the cross sections
    
    # Translation of flow energy to speed for future use:
    
    import numpy as np
    V_array = np.sqrt(2*Energy_array_J/m)      # [m/s] Speeds associated with the cross sections
     
    # Generation of the vector of velocities where the Maxwellian is evaluated:
    
    v = np.linspace(0, 1*V_array.max(), 200)  # [m/s] vector of speeds to increase the accuracy
    
    # Initialization of the generated arrays:
    
    Maxwellian_distribution = np.zeros(len(v))
    product = np.zeros(len(v))
    
    # Integration process:
    
    for i in range(len(v)):
        Maxwellian_distribution[i] = 4*np.pi*((m/(2*np.pi*T))**(3/2))*(v[i]**2)*np.exp(-m*(v[i]**2)/(2*T))   # [s/m] Maxwellian differential density
        # As the cross section data is only provided inside a small range, it must be interpolated outside it:
        if v[i]<V_array.min():
            local_cross = 0 # below the threshold energy collisions dont take place
        elif v[i]>=V_array.min() and v[i]<=V_array.max():
            local_cross = np.interp(v[i], V_array, CS_array) # Inside the region, the cross section is interpolated
        elif v[i]>V_array.max():
            f_cs = CS_array[-2] + (CS_array[-1]-CS_array[-2])*(v[i] - V_array[-2])/(V_array[-1]-V_array[-2]) # outside the range, the cross section is linearly extrapolated
            if f_cs<=0:
                local_cross=0 # it cannot be negtive
            else:
                local_cross=f_cs    
        product[i] = local_cross*v[i]*Maxwellian_distribution[i] # [m^2] diferential collision rate coefficient
        
    product_integr = np.trapz(product, v)       # [m^3/s] integral over the velocity range
    n = np.trapz(Maxwellian_distribution, v)    # the density is included because at it is only evaluated at some parts it is < 1 and 
                                                # compensates for the regions neglected
    
    collision_rate = product_integr/n           # [m^3/s] collision rate coefficient
    
    return collision_rate
    
#############################################################################################################
#############################################################################################################
    
def collision_rate_BiMaxwellian(Database_root, Type, specie1, specie2, process, file_name, Vel_par, Vel_perp):
    
    """
    This codes computes de collision rate assuming an an-isotropic system with  
    a given electron flow velocity in the parallel and perpendicular directions.
            
    INPUTS REQUIRED:
    
    -Database_root: the root directory where the database is stored (in the
                    desktop, in a specific folder, in a USB...)
    -Type: type of the desired interaction. With the current state of the database
          this field must be one of: 'BINARY' or 'UNARY'
    -specie1: first of the species involved in the interaction
    -specie2: second specie involved in the interaction
        
        NOTE_1: the species are introduced with their sign an their electronic state.
                Some examples are: 'e-' for electron
                                   'Xe' for neutral Xe
                                   'Xe+2' for doble ionized Xenon--> Symbol+Charge+#charges
                                   'Xe[5s2]' for excited Xe --> Symbol[electronic state of the las valence e-]
                                   
        NOTE_2: the specied dont necessarily need to be entered in order.
        
        NOTE_3: If the type UNARY is introduced, the secondary specie should be left as ''
                                         
    -process: desired reaction. 'i1' for one ionization of 1 e-
                                'elastic' for elastic collision
                                'e[...]' for excitation to the [...] electronic configuration
                                
    -file_name: name of the file that wants to be read (it must be a .txt)
    
        NOTE_3: If the name of the .txt is unknown, leave it as '' and the code
                will automatically read the first file associated with a database or,
                in the case where only one file is stored, it will take read this file.
                
    -Vel_par: velocity of the flow of electrons in the parallel direction. In [eV]
    -Vel_perp: velocity of the flow of electrons in the perpendicular direction. In [eV]

    OUTPUT:
    
    -The collision rate at the specified energies. In [m^3/s]
    """
    
    # Dictionary with collision data
    collision_info = get_collision_info(Database_root, Type, specie1, specie2, process, file_name)
    
    # Extraction of energy and differential cross section values from dictionary

    [Energy_array, CS_array] =_get_E_and_CS_arrays_from_info(collision_info) #[eV, m^2]
    
    # Important parameteres for the Maxwellian distribution:
    
    m = 9.10938356e-31              # [Kg] electron mass
    K = 1.6e-19                     # [J/eV] Boltzmann constant
    
     # Adjusting of the flow velocity to S.I.:
    
    T_par = K*Vel_par               # [J] Flow temperature in the parallel direction
    T_perp = K*Vel_perp             # [J] Flow temperature in the perpendicular direction
    Energy_array_J = K*Energy_array       # [J] Energies associated with the cross sections
     
    # Translation of flow energy to speed for future use:
    
    import numpy as np
    V_array = np.sqrt(2*Energy_array_J/m)      # [m/s] Speeds associated with the cross sections
    
   # Generation of the vectors of velocities where the Bi-Maxwellian is evaluated:

    V_par = np.linspace(0, V_array.max(), 200) # [m/s] Vector of velocities in the parallel direction
    V_perp = np.linspace(0, V_array.max(), 200)# [m/s] Vector of velocities in the perpendicular direction 
    
   # Initialization of the generated arrays:

    BiMaxwellian_matrix_prob = np.zeros((len(V_perp),len(V_perp)))
    matrix_product = np.zeros((len(V_par),len(V_perp)))
    
    # Integration process:

    for pos_row in range(len(V_par)):
        v_par = V_par[pos_row]
        for pos_col in range(len(V_perp)):
            v_perp = V_perp[pos_col]
            BiMaxwellian_matrix_prob[pos_row, pos_col] = ((m**(3/2))/(np.sqrt(2*np.pi*T_par)*T_perp))*v_perp*np.exp(-m*(v_par**2)/(2*T_par) - m*(v_perp**2)/(2*T_perp))   # [s/m]            
            V = np.sqrt(v_par**2 + v_perp**2)
            # As the cross section data is only provided inside a small range, it must be interpolated outside it:
            if V<V_array.min():
                local_cross = 0 # below the threshold energy collisions dont take place
            elif V>=V_array.min() and V<=V_array.max():
                local_cross = np.interp(V, V_array, CS_array) # Inside the region, the cross section is interpolated
            elif V>V_array.max():
                f_cs = CS_array[-2] + (CS_array[-1]-CS_array[-2])*(V - V_array[-2])/(V_array[-1]-V_array[-2]) # outside the range, the cross section is linearly extrapolated
                if f_cs<=0:
                    local_cross=0 # it cannot be negtive
                else:
                    local_cross=f_cs  
            matrix_product[pos_row, pos_col] = local_cross*V*BiMaxwellian_matrix_prob[pos_row, pos_col] # [m^2] diferential collision rate coefficient
  
    # integration over the parallel flow direction:
    integration_parallel = 2*np.trapz(matrix_product, V_par)
    # integration over the perpendicular flow direction
    integration_par_perp = np.trapz(integration_parallel, V_perp) # [m^3/s] integral over the velocity range
      
    density_parallel = 2*np.trapz(BiMaxwellian_matrix_prob, V_par)
    n = np.trapz(density_parallel, V_perp) # the density is included because at it is only evaluated at some parts it is < 1 and 
                                           # compensates for the regions neglected
    
    collision_rate = integration_par_perp/n # [m^3/s] collision rate coefficient
        
    return collision_rate

#############################################################################################################
#############################################################################################################
    
    
    
def _select_fileroot_filename_unknown(root,Type, Specie1, Specie2, reaction): 
    
    [Database, tree] = _get_directory_structure(root)
    if Type == 'BINARY':
        if len(Database['DATABASE'][Type][Specie1][Specie2][reaction]) == 1:
            for file_name in Database['DATABASE'][Type][Specie1][Specie2][reaction].keys():#database.txt,...
                file_root = root+'/'+Type +'/'+Specie1+'/'+Specie2+'/'+reaction+'/'+file_name
        else:
            options = []
            for file_name in Database['DATABASE'][Type][Specie1][Specie2][reaction].keys():#database.txt,...
                if len(Database['DATABASE'][Type][Specie1][Specie2][reaction]) != 1 and file_name[0:8] == 'database':
                    options = options + [file_name]
            file_used = options[0]
            file_root = root+'/'+Type +'/'+Specie1+'/'+Specie2+'/'+reaction+'/'+file_used

    elif Type == 'UNARY':
        for specie1 in Database['DATABASE'][Type].keys():# e-, Ar, Kr,...
            if specie1==Specie1:
                for reaction in Database['DATABASE'][Type][specie1].keys():#de_ex
                    if len(Database['DATABASE'][Type][specie1][reaction]) == 1:
                        for file_name in Database['DATABASE'][Type][specie1][reaction].keys():#database.txt,...
                            file_root = root+'/'+Type +'/'+specie1+'/'+reaction+'/'+file_name
                    else:
                        options = []
                        for file_name in Database['DATABASE'][Type][specie1][reaction].keys():#database.txt,...
                            if len(Database['DATABASE'][Type][specie1][reaction]) != 1 and file_name[0:8] == 'database':
                                options = options + [file_name]
                        file_used = options[0]
                        file_root = root+'/'+Type +'/'+specie1+'/'+reaction+'/'+file_used
                            
        
    return file_root

###############################################################################################################
###############################################################################################################

def _get_directory_structure(rootdir):
    
    """
    This function reads the entire database and strutures it as a DICTIONARY and a TREE
    for the ease of looking to what the Database currently stores:
        
        INPUTS: rootdir = root directory where the database is stored in the computer, USB,... Ex: rootdir = 'C:/Users/Antonio/Desktop/DATABASE_YAML'
        
        OUTPUTS: Dictionay = whole database organized as a dictionary variable-->there are files, subfiles, subsubfiles...
                 Tree = whole database printed as a tree --> it can be easilly checked the current information stored in the database
    """
    
    import os
    import yaml
    from functools import reduce

    # First the disctionary structure where the database will be stored is created:
    
    dir = {}
    
    # Then the directory provided is covered and read:
    
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        subdir = dict.fromkeys(files)
        parent = reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
    
    # The output variables as dictionary and tree are specified    
    
    Dictionary = dir
    Tree = yaml.dump(Dictionary, default_flow_style = False)
    
    return [Dictionary, Tree]

###############################################################################################################
###############################################################################################################
def _get_E_and_CS_arrays_from_info(Information):
    
    """
    This function obtains from the information provided by the user 
    the Cross_section values IN S.I. (m^2) associated to different Energy values IN eV
    
    INPUTS: Information = corresponds to the information stored in a file of the database which has previously been extracted with "fun_get_collision_info"
    
    OUTPUTS: Energy_array = [eV], array of energies for which the cross section value is provided
             Cross_section_array = [m^2] array of the cross sections associated to each of the energies. 
    """
    # Depending of the database or model used to obtain the data, the cross sections are obtained in a different way:
    
    if Information['TYPE'] == 'Database': # The data was obtained from a database
        Energy_array, Cross_section_array = _from_database(Information['DATA'])
    elif Information['TYPE'] == "Model-Drawin" :  # The Drawin Model was followed
         Energy_array, Cross_section_array = _from_Drawin(Information['DATA'])
    elif Information['TYPE'] == 'Model-Gryzinski' : # Grazynski Model was follwed
         Energy_array, Cross_section_array = _from_Gryzinski(Information['DATA'])
        
    return Energy_array, Cross_section_array

def _from_database(DATA):
    
    """
    This function extracts the cross section and energy values from the files 
    that store information obtained from a database
    """
    import numpy as np
    
    energy = np.array(DATA['ENERGY']['VALUES'])
    cross_section = np.array(DATA['CROSS-SECTION']['VALUES'])
    if DATA['CROSS-SECTION']['UNITS'] == 'cm^2':
        cross_section = cross_section*1e-4
    return(energy, cross_section)  

def _from_Drawin(DATA):
    
    """
    This function computes the cross section and energy values from the files 
    that store information following the Drawin Model
    """
    import numpy as np
    
    a_0 = DATA['a_0']['VALUES']
    epsilon_i_H = DATA['epsilon_i_H']['VALUES']
    epsilon_i = DATA['epsilon_i']['VALUES']
    xi = DATA['xi']['VALUES']
    beta_1 = DATA['beta_1']['VALUES']
    beta_2  =DATA['beta_2']['VALUES']
    final_E = DATA['Final_E']['VALUES']    
    Energy_range = np.linspace(epsilon_i, final_E, 200)
     
    u = Energy_range/epsilon_i
        
    g = ((u-1)/(u**2))*np.log(1.25*beta_2*u)
        
    Cross_sections = 2.66*np.pi*(a_0**2)*beta_1*((epsilon_i_H/epsilon_i)**2)*xi*g

    return(Energy_range, Cross_sections)

def _from_Gryzinski(DATA):
    
    """
    This function computes the cross section and energy values from the files 
    that store information following the Gryzinski Model
    """
    import numpy as np
    
    a_0 = DATA['a_0']['VALUES']
    epsilon_i_H = DATA['epsilon_i_H']['VALUES']
    epsilon_i = DATA['epsilon_i']['VALUES']
    xi = DATA['xi']['VALUES']
    final_E = DATA['Final_E']['VALUES']    
    Energy_range = np.linspace(epsilon_i, final_E, 200)
     
    u = Energy_range/epsilon_i
    
    gg = (1+2/3*(1-1/(2*u))*np.log(np.e+(u-1)**(1/2)))    
    g = ((u-1)/u**2)*((u/(u+1))**(3/2))*((1-1/u)**(1/2))*gg
        
    Cross_sections = 4*np.pi*(a_0**2)*((epsilon_i_H/epsilon_i)**2)*xi*g

    return(Energy_range, Cross_sections)
        
###############################################################################################################
###############################################################################################################

