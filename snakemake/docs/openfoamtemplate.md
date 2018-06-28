[Home](../README.md) 

# Setting up a new OpenFOAM #

## Create a new folder on the fileshare ##

Create a new folder on the fileshare for your OpenFOAM analysis and create a folder for the Batch-Shipyard configuration files:

~~~~
cd $FILESHARE
mkdir <yourfolder>
cd <yourfolder>
mkdir azurebatch
cd azurebatch
~~~~

Batch-shipyard needs 4 configuration files to be able to create pools and run jobs. These can be copied from the ones below. Batch-shipyard also has comprehensive documentation: https://github.com/Azure/batch-shipyard/tree/master/docs. There are also examples of the full settings available: https://github.com/Azure/batch-shipyard/tree/master/config_templates


## config.yaml ##

~~~~
batch_shipyard:
  storage_account_settings: mystorageaccount
global_resources:
  docker_images:
  - alfpark/openfoam:4.0-gcc-openmpi
  volumes:
    shared_data_volumes:
      sharedfiles:
        volume_driver: azurefile
        storage_account_settings: mystorageaccount
        azure_file_share_name: openfoam
        container_path: $AZ_BATCH_NODE_SHARED_DIR/gfs
        mount_options:
        - file_mode=0777
        - dir_mode=0777
~~~~

If using a fileserver modify the shared shared_data_volumes to point to the storage_cluster:

~~~~
  volumes:
    shared_data_volumes:
      mystoragecluster:
        volume_driver: storage_cluster
        container_path: $AZ_BATCH_NODE_SHARED_DIR/mystoragecluster
        mount_options: []
~~~~

## credentials.yaml ##

~~~~
credentials:
  batch:
    account_key: youraccountkeyendingin==
    account_service_url: https://hpcuoaopenfoambatch.australiasoutheast.batch.azure.com/
  storage:
    mystorageaccount:
      account: hpcuoaopenfoam
      account_key: Oyourstorageaccountkeyendingin==
      endpoint: core.windows.net
~~~~

## jobs.yaml ##

~~~~
job_specifications:
  - id: openfoam
    auto_complete: true
    tasks:
    - docker_image: alfpark/openfoam:4.0-gcc-openmpi
      shared_data_volumes:
      - sharedfiles
      multi_instance:
        num_instances: pool_current_dedicated
      command: $AZ_BATCH_NODE_SHARED_DIR/gfs/run.sh
~~~~


## pool.yaml ##

The vm_size can be modified here. For a full list go to: https://docs.microsoft.com/en-us/azure/virtual-machines/windows/sizes

~~~~
pool_specification:
  id: docker-openfoam-tcp1
  vm_configuration:
    platform_image:
      offer: UbuntuServer
      publisher: Canonical
      sku: 16.04-LTS
  vm_count:
    dedicated: 2
    low_priority: 0
  vm_size: Standard_D1
  inter_node_communication_enabled: true
  ssh:
    username: hpcadmin
~~~~

## Setting up a fileserver cluster ###

To setup a fileserver to use instead of an Azure Blob storage account an additional configuration file is needed 'fs.yaml'. This config defines the specifications of the file server eg: number of servers/disk etc. Using the fileserver for Azure Batch jobs is very similar to using a storage account. 

https://github.com/Azure/batch-shipyard/blob/master/docs/65-batch-shipyard-remote-fs.md

*fs.yaml*

~~~~
remote_fs:
  resource_group: hpcuaopenfoam
  location: australiasoutheast
  managed_disks:
    premium: false
    disk_size_gb: 128
    disk_names:
    - p30-disk0a
    - p30-disk1a
  storage_clusters:
    mystoragecluster:
      hostname_prefix: hpcuaopenfoamsc
      ssh:
        username: shipyard
      file_server:
        mount_options:
        - noatime
        - nodiratime
        mountpoint: /data
        type: nfs
      network_security:
        ssh:
        - '*'
      virtual_network:
        address_space: 10.0.0.0/16
        existing_ok: true
        name: hpcuaopenfoamscvnet
        subnet:
          address_prefix: 10.0.0.0/24
          name: hpcuaopenfoamsc-server-subnet
      public_ip:
        enabled: true
        static: false
      vm_count: 1
      vm_size: Standard_F1s
      vm_disk_map:
        '0':
          disk_array:
          - p30-disk0a
          - p30-disk1a
          filesystem: btrfs
          raid_level: 0
          
~~~~

Setting up the fileshare cluster:

~~~~
SHIPYARD_CONFIGDIR=. $SHIPYARD/shipyard fs disks add
~~~~

This will SSH into the first (and only) VM in the storage cluster:

~~~~
SHIPYARD_CONFIGDIR=. $SHIPYARD/shipyard fs cluster add mystoragecluster
~~~~

To delete the file server, you can issue:

~~~~
 # keep the data disks, resource group and virtual network
SHIPYARD_CONFIGDIR=. $SHIPYARD/shipyard fs cluster del mystoragecluster
~~~~
	

	 
		 
		 
