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
########### INPUTS REQUIRED ###########################

Database_root = 'C:/Users/Antonio/Desktop/DATABASE'
Type = 'BINARY'
specie1 = 'e-'
specie2 = 'Xe'
process = 'i1'
file_name = ''

########### CALL TO THE DATABASE  #####################

import db_codes
Collision_info = db_codes.get_collision_info(Database_root, Type, specie1, specie2, process, file_name)
