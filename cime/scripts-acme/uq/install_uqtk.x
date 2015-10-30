#!/bin/bash



## Open the UQ-workflow tarball                     
#tar -xvzf uq_workflow.tar                                      

## Open the UQTk tarball                             
tar -xvzf uqtk_v2.2.tgz                                        

## Config UQTk by creating config.<machinename> file 
### and linking to or copying it to config.site                 
### see the example config files                      
cd UQTk_v2.2/config                     
ln -sf config.gnu config.site                                  

## Make and cd up to the main folder                 
cd ..
make
cd ..

