# Cross-section-database-for-hybrid-PIC-fluid-plasma-codes
This database was created with the purpose of providing easy acces to relevant collision cross sections and rate coefficients for hybrid PIC/fluid plasma codes
This repository contains: 

- Licence.

- folder DATABASE: which contains the folder tree stucture with all the files of the current collision cross sections stored in the database. Inside this database the are two sections. To navigate inside the database, the most electronegative particle is identified in the first place, then the colliding partner, type of collision process and the .txt file containing the cross section values. The type of processes currently stored in the database are:

          - i1: ionization of one electron (one electron is freed).
          - i2: ionization of two electron.
          - elastic: elastic collision between particles.
          - e[]: excitation reaction....e[electronic configuration of the last valence electron of the product]
          - ch_e: charga exchange interaction
          
   In addition, the name of the .txt file indicate the source and the name of the source. This information is also found inside each of the files, with a detailed explanation of the reaction.
       
- A .xlsx file showing the current status of the database in a grafical way. This provides a direct view of what collisions are alrready stored.

- folder DATABASE_CODES: which contains the complementary codes of the database to allow communicotion with hybrid PIC/fluid codes. Inside this folder all the code functions are stored in the file: db_codes.py . Four main functions are designed to be called by the user:                 - get_collision_info
        - interpolate_at_velocity
        - collision_rate_Maxwellian
        - collision_rate_BiMaxwellian
        
   Additionally, there are four sample files that describe the four main codes designed to interact with the user of the hybrid code. These describe the objective of the function and provide an explanation of the required inputs and the provided outputs.
