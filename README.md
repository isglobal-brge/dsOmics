# Package Description
This package aims to implement omic data anlalysis using DataSHIELD infrastructure

# Installation
The following steps describe how to install the DataSHIELD training environment to test the first version of the function for methylation analysis:


1. Here are the instructions on how to install the virtual box and the two training Virtual Machines:
https://data2knowledge.atlassian.net/wiki/spaces/DSDEV/pages/12943427/Instructions+for+Windows+and+Mac+users#InstructionsforWindowsandMacusers-DownloadandinstallVirtualBox{:target="_blank"}

2. Here are the instructions on how to create a new project in your Virtual machines. (For simplicity, I have called this project ‘OMICS' in my virtual machines):
https://data2knowledge.atlassian.net/wiki/spaces/DSDEV/pages/12943489/Using+your+own+data+in+Opal#UsingyourowndatainOpal-CreatinganewOpalproject

You can use the attached .csv and .xlsx files that are included in the folders server1 and server2, to upload the data and their dictionaries to the new project in the two opal servers.

3. Then, you have to add the server-side function (lmFeatureDS) in both servers following these instructions:
https://data2knowledge.atlassian.net/wiki/spaces/DSDEV/pages/12943483/How+to+upload+a+server-side+function


4. You need also to install the dsBetaTest package using the following instructions:
https://data2knowledge.atlassian.net/wiki/spaces/DSDEV/pages/12943493/Update+or+install+packages+on+the+training+Opal+servers
and to add the new disclosure-control filters by login to each server and under Administration/DataSHIELD/Options 


