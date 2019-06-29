"""
CODE DESCRIPTION: 
    
    This code exemplifies how the database would be called to extract the file
    with all the information of an specified collision.
    
INPUTS REQUIRED:
    
    -Database_root: the root directory where the database is stored (in the
                    desktop, in a specific folder, in a USB...)
    -Type: type of the desired interaction. With the current state of the database
          this field must be one of: 'BINARY' or 'UNARY'
    -specie1: first of the species involved in the interaction
    -specie2: secod specie involved in the interaction
        
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
                
    -Velocity: velocity at which the differential cross-section is desired. In [m/s]

OUTPUT:
    
    -The cross section of the specified reation at the specified velocity. In [m^2]
"""
########### INPUTS REQUIRED ###########################

Database_root = 'C:/Users/Antonio/Desktop/DATABASE'
Type = 'BINARY'
specie1 = 'Xe'
specie2 = 'Xe+1'
process = 'ch_e'
file_name = ''

Velocity = 383.13 # [m/s]

########### CALL TO THE DATABASE  #####################

import db_codes
interpolated_differential_cross_section = db_codes.interpolate_at_velocity(Database_root, Type, specie1, specie2, process, file_name, Velocity) #[m^2]
