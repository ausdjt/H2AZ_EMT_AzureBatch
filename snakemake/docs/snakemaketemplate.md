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
    - hpcanudocker.azurecr.io
  docker_images:
  - hpcanudocker.azurecr.io/anu:latest
  volumes:
    shared_data_volumes:
      mystoragecluster:
        volume_driver: storage_cluster
        container_path: /data
        mount_options: []
        bind_options: rw
~~~~

## credentials.yaml ##

~~~~
credentials:
  batch:
    aad:
      endpoint: https://batch.core.windows.net/
      directory_id: aaddirectoryguidfromtheportal
      user: username@yourtennant
      password: abcpassword
      token_cache:
        enabled: true
        filename: .aad_token_cache1
    account_service_url: https://yourbatchaccount.australiasoutheast.batch.azure.com/
    resource_group: yourresourcegroupname
  management:
    aad:
      endpoint: https://management.azure.com/
      directory_id: aaddirectoryguidfromtheportal
      user: username@yourtennant
      password: abcpassword
      token_cache:
        enabled: true
        filename: .aad_token_cache
    subscription_id: yousubscriptionguidfromtheportal
  storage:
    mystorageaccount:
      account: yourstorageaccountname
      account_key: yourstorageaccountkeyendingin==
      endpoint: core.windows.net
  docker_registry:
    youdockerregistry.azurecr.io:
      username: dockeruser
      password: dockerpasswordfromtheportal
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
  - id: snakemake
    auto_complete: true
    user_identity:
      specific_user:
        gid: 1001
        uid: 1001
    tasks:
    - docker_image: hpcanudocker.azurecr.io/anu:latest
      shared_data_volumes:
      - mystoragecluster
      command: /data/jobrun.sh
~~~~


## pool.yaml ##

The vm_size can be modified here. For a full list go to: https://docs.microsoft.com/en-us/azure/virtual-machines/windows/sizes

~~~~
pool_specification:
  id: snakemake
  virtual_network:
    name: hpcanusnakemakescvnet
    resource_group: hpcanusnakemake
    address_space: 10.0.0.0/16
    subnet:
      name: hpcanusnakemakesc-server-subnet
      address_prefix: 10.0.0.0/24
  vm_configuration:
    platform_image:
      offer: UbuntuServer
      publisher: Canonical
      sku: 16.04-LTS
  vm_count:
    dedicated: 1
    low_priority: 0
  vm_size: Standard_E2_v3
  inter_node_communication_enabled: false
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
		 
		 
		 
