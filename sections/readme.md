# Section Generation
Create section files contained the QCed data in the google drive processed directory.  
All the "heavy lifting" is done by the `Instrument` subclasses which contain all the properties and methods required to extract the data from the processed directory into a matlab structure.

## Scripts
- *run_generate_sections* shell script that loads the matlab module and runs *generate_sections.m*
- *generate_sections.m* matlab script that creates sections as defined in the *processed/survey_metadata* directory of the google drive.

## Classes
- *Vessel.m* matlab class that handles all the instruments for a vessel.  
  **Properties**
  - `instruments`: cell array of `Instrument` subclasses

  **Methods**
  - `get_all_data(start,stop)`:
      calls the `get_data` method of every instrument and returns a structure with data for that vessel between `start` and `stop`  
 
## Instruments
Subdirectory containing the `Instrument` class and subclasses for individual implementations for different instruments
- *Instrument.m* abstract class that defines the protocol for individual instruments to follow  
  **Properties**
  - `name`: `string` that defines the instrument name. Is used to define the variable name in the section files.
  - `variables`: `string array` of variables names to include in the section files. If, as by default, this is the empty string array (not the empty string)
                  then every variable should be included.
  - `data_source` (abstract): a file or directory path to the data for the particular instrument

  **Methods (Abstract)**
  - `get_data(start,stop)`: returns a structure or cell array of structures with fields defined by `variables` with the data between `start` and `stop`  

- *ADCP_Ship_Combo* implementation (subclass) of `Instrument` for the shipboard ADCP combo files.   
  **Constructor Arguments**  
  - `name`: name of this ADCP in the section files e.g. *ADCP_PE_wh1200*
  - `datafile`: path to the file
  - `variables`: `string array` of variables to include. Default: `string.empty` in which case all variables are included.

- *Hydro_Combo* implementation of `Instrument` for the combined VMP and CTD files.  
  **Constructor Arguments**  
  - `name`: name in the section files e.g. *HYDRO_Pelican*
  - `data_file`: path to the hydro combo file e.g. *'.../Processed/HydroCombo/SUNRISE2021_PE_HydroCombo_Processed.mat'*
  - `variables`: `string array` of variables to include. Default: `string.empty` in which case all variables are included.



- *TChain* implementation of `Instrument` for the TChains.  
  **Constructor Arguments**
  - `name`: name in the section files e.g. *TCHAIN_Pelican*, n.b. if the section contains more than one TChain deployment then `Vessel.get_all_data` will turn the cell array of structures returned from `get_data` into multiple structures with names *name_1*, *name_2*, ...
  - `data_directory': path to the directory containing all the TChain deployments for a single vessel
  - `variables`: `string array` of variables to include. Default: `string.empty` in which case all variables are included.

- *ADCP_Rhib* redundant implementation of `Instrument` for the rhib adcps.  
  This class used the intermediate files rather than the processed files and will need reimplementing once those files are uploaded to the google drive.

 - *VMP_Combo* redundant implementation of `Instrument` for the VMP combo files.  
  VMP data is merged with CTD into Hydro files
