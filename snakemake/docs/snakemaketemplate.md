[Home](../README.md) 

# Setting up a new Snakemake #

## Create a new folder on the fileshare ##

Create a new folder on the fileshare for your Snakemake workflow and create a folder for the Batch-Shipyard configuration files:

~~~~
cd $FILESHARE
mkdir <yourfolder>
cd <yourfolder>
mkdir azurebatch
cd azurebatch
~~~~

Batch-shipyard needs 4 configuration files to be able to create pools and run jobs. These can be copied from the ones below. Batch-shipyard also has comprehensive documentation: https://github.com/Azure/batch-shipyard/tree/master/docs. There are also examples of the full settings available: https://github.com/Azure/batch-shipyard/tree/master/config_templates.

These template files are also available in the GitHub repository.


## config.yaml ##

~~~~
batch_shipyard:
  storage_account_settings: mystorageaccount
global_resources:
  additional_registries:
    docker:
    - hpcuoadocker.azurecr.io
  docker_images:
  - hpcuoadocker.azurecr.io/rnaseq:latest
  volumes:
    shared_data_volumes:
      sharedfiles:
        volume_driver: azurefile
        storage_account_settings: mystorageaccount
        azure_file_share_name: fileshare
        container_path: /home/hpcadmin/fileshare
        mount_options:
        - file_mode=0777
        - dir_mode=0777


~~~~

## credentials.yaml ##

~~~~
credentials:
  batch:
    account_key: yourbatchaccountkeyendingin==
    account_service_url: https://hpcuoasnakemakebatch.australiasoutheast.batch.azure.com/
  storage:
    mystorageaccount:
      account: hpcuoasnakemake
      account_key: yourstorageaccountkeyendingin==
      endpoint: core.windows.net
 ~~~~

To use an Azure Container registry add the following configuration: 

~~~~
 docker_registry:
    hpcuoadocker.azurecr.io:
      username: hpcuoadocker
      password: yourdockerlogin
~~~~
## jobs.yaml ##

~~~~
job_specifications:
  - id: rnaseqjoball
    tasks:
    - docker_image: hpcuoadocker.azurecr.io/rnaseq:latest
      shared_data_volumes:
      - sharedfiles
      command: /home/hpcadmin/fileshare/RNAseq_snakemake-masterAll/jobrunall.sh
~~~~


## pool.yaml ##

The vm_size can be modified here. For a full list go to: https://docs.microsoft.com/en-us/azure/virtual-machines/windows/sizes

~~~~
pool_specification:
  id: rnaseqpoolall
  vm_configuration:
    platform_image:
      offer: UbuntuServer
      publisher: Canonical
      sku: 16.04-LTS
  vm_count:
    dedicated: 0
    low_priority: 1
  vm_size: Standard_D1
  ssh:
    username: hpcadmin
~~~~

## Snakemake Shell command ##

For the Snakemake steps use this basic template:

~~~~
         #!/usr/bin/env bash
cd /data
#snakemake --latency-wait 60 --snakefile ./Development/H2AZ_EMT/snakemake/workflows/MDCK_ChIP-Seq.py --configfile ./Development/H2AZ_EMT/snakemake/configs/config.json --config ASSAY=ChIP-Seq RUNID=NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq WORKFLOWDIR=Development --jobs -pr
snakemake --latency-wait 60 --snakefile ./Development/H2AZ_EMT/snakemake/workflows/MDCK_RNA-Seq.py --configfile ./Development/H2AZ_EMT/snakemake/configs/config.json --config ASSAY=RNA-Seq RUNID=NB501086_0082_RDomaschenz_JCSMR_mRNAseq WORKFLOWDIR=Development --jobs -pr

~~~~
		 
		 
		 
