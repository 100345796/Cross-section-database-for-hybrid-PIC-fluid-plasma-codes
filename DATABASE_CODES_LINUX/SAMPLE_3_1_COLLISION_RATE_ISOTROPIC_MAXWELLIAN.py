"""
CODE DESCRIPTION: 
    
    This code exemplifies how the database would be called to compute the collision
    rate for a specific reaction at a given temperature
    
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
########### INPUTS REQUIRED ###########################

Database_root = 'C:/Users/Antonio/Desktop/DATABASE'
Type = 'BINARY'
specie1 = 'e-'
specie2 = 'Kr'
process = 'i1'
file_name = ''

Velocity = 20 # [eV]

########### CALL TO THE DATABASE  #####################

import db_codes
Collision_rate = db_codes.collision_rate_Maxwellian(Database_root, Type, specie1, specie2, process, file_name, Velocity) #[m^3/s]


