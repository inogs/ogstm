
         trn(:,:,:,ppO2o) = 250.0

* pelagic nutrients (mMol /m3)
         trn(:,:,:,ppN1p) = 0.075
         trn(:,:,:,ppN3n) = 3.0
         trn(:,:,:,ppN4n) = 0.5
         trn(:,:,:,ppO4n) = 0.5
         trn(:,:,:,ppN5s) = 5.5
         trn(:,:,:,ppN6r) = 1.

* pelagic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppR6c) = 17.0
         trn(:,:,:,ppR6n) = 0.24
         trn(:,:,:,ppR6p) = 0.02
         trn(:,:,:,ppR6s) = 0.10

* dissolved organic matter

         trn(:,:,:,ppR1c) = 0.0001
         trn(:,:,:,ppR1n) = trn(:,:,:,ppR1c)/12./6.7
         trn(:,:,:,ppR1p) = trn(:,:,:,ppR1c)/12./106.
CCC F79 20 10 20004 we put trn(:,:,:,ppR1s) = 0.1 
CCC      trn(:,:,:,ppR1s) = 0.
         trn(:,:,:,ppR1s) = trn(:,:,:,ppR1c)/12./6.7
         trn(:,:,:,ppR2c) = trn(:,:,:,ppr1c)
         trn(:,:,:,ppR7c) = trn(:,:,:,ppr1c)

*-----------------------------------------------------
* State variables for phytoplankton model
*-----------------------------------------------------

* pelagic diatoms  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppP1c) = 8.
         trn(:,:,:,ppP1n) = 0.114
         trn(:,:,:,ppP1p) = 0.018
         trn(:,:,:,ppP1s) = 0.026
         trn(:,:,:,ppP1i) = 0.1

* pelagic flagellates  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppP2c) = 5.9
         trn(:,:,:,ppP2n) = 0.0926
         trn(:,:,:,ppP2p) = 0.0076
         trn(:,:,:,ppP2i) = 0.1

* picophytoplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppP3c) = 5.9
         trn(:,:,:,ppP3n) = 0.0926
         trn(:,:,:,ppP3p) = 0.0076
         trn(:,:,:,ppP3i) = 0.1

* pelagic inedibles  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppP4c) = 5.9
         trn(:,:,:,ppP4n) = 0.0926
         trn(:,:,:,ppP4p) = 0.0076
         trn(:,:,:,ppP4i) = 0.1

*-----------------------------------------------------
* State variables for mesozooplankton model
*-----------------------------------------------------
* mesozooplankton have a fixed ratio n and p in this model:

* carnivorous mesozooplankton ( mg C/m3 )

         trn(:,:,:,ppZ3c) = 1.2
         trn(:,:,:,ppz3n) = trn(:,:,:,ppz3c)/12./6.7
         trn(:,:,:,ppz3p) = trn(:,:,:,ppz3c)/12./106.

* omnivorous mesozooplankton ( mg C/m3 )

         trn(:,:,:,ppZ4c) = 1.2
         trn(:,:,:,ppz4n) = trn(:,:,:,ppz4c)/12./6.7
         trn(:,:,:,ppz4p) = trn(:,:,:,ppz4c)/12./106.

*-----------------------------------------------------
* State variables for microzooplankton model
*-----------------------------------------------------

* pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppZ5c) = 7.2
         trn(:,:,:,ppZ5n) = 0.12
         trn(:,:,:,ppZ5p) = 0.0133

* heterotrophic flagellates (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppZ6c) = 2.421
         trn(:,:,:,ppZ6n) = 0.0508
         trn(:,:,:,ppZ6p) = 0.00470

*-----------------------------------------------------
* State variables for pelagic bacteria model B1
*-----------------------------------------------------
* pelagic bacteria  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

         trn(:,:,:,ppB1c) = 15.7
         trn(:,:,:,ppB1n) = 0.26
         trn(:,:,:,ppB1p) = 0.029

