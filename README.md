# SurfFlow<a name="toc"></a>

SurfFlow is a collection of workflows designed to generate surfaces with arbitrary Miller indices and terminations and perform VASP calculations to
automatically calculate surface energies and Wulff shapes. surfflow uses the [FireWorks](https://materialsproject.github.io/fireworks/index.html)
framework and relies heavily
on [pymatgen](https://pymatgen.org/), [atomate](https://atomate.org/), and the [Materialsproject](https://materialsproject.org/) in general.

1. [Download and Installation](#install)
    1. [Create a virtual environment for surfflow](#virtualenv)
    2. [Install and configure MongoDB locally if you want your database to run also in conda](#mongodb)
    3. [Create a folder structure, download and install](#surfflowinstall)
    4. [Configuring FireWorks](#configurefw)
    5. [Configuring pymatgen](#configurepymatgen)
2. [Testing your installation of FireWorks](#testing)
3. [A note about atomate](#customatomate)
4. [Using surfflow](#using)
    1. [Running a workflow](#runningawf)
    2. [Looking at the results](#results)

## Download and Installation <a name="install"></a>

We recommend installing surfflow inside a virtual python environment. We use [mamba](https://github.com/mamba-org/mamba) in this guide to create and manage environments,
but [conda](https://docs.conda.io/en/latest/), [virtualenv](https://virtualenv.pypa.io/en/latest/), or similar are also possible. We currently
recommend **mamba** because it is essentially **conda**, but with much faster environment solving and package installation. If you prefer to proceed with
conda instead, you can just replace all the `mamba` commands with `conda`. In fact, in the later sections of this guide, we use mamba and conda
interchangeably in some commands and when referring to environments.

### Create a virtual environment for SurfFlow <a name="virtualenv"></a>

1. Make sure that you have mamba installed or download and install mambaforge from [here](https://github.com/conda-forge/miniforge). If you already
   have conda installed, you can also install mamba by running `conda install mamba -n base -c conda-forge`.
2. Update mamba and then create a new python 3 environment named "surfflow" (or whatever you decide) by typing `mamba update mamba` followed
   by `mamba create --name surfflow python=3`
3. Switch to your new environment: `mamba activate surfflow`

[Back to top](#toc)

### Install SurfFlow from PyPi <a name="pypi"></a>

Simply install surfflow from PyPi with the following command

`pip install surfflow`

**Note**: `atomate` is a dependency of the `surfflow` package. However, due to some changes that either fix certain bugs or add necessary functionality yet to be implemented in the official atomate version, we use a slightly modified fork of `atomate`. This fork should be installed with

`pip install git+atomate@git+https://github.com/fyalcin/atomate.git`

### Install and configure MongoDB locally if you want your database to run also in conda environment<a name="mongodb"></a>

This is **optional**, and you should skip these steps if you **already have a database running** somewhere that you want to use and can access from
the machine you are using. See also [this section](https://atomate.org/installation.html#mongodb) of the atomate installation instructions.

1. Install the mamba package of MongoDB by simply running `mamba install mongodb` while the `surfflow` environment is active.
2. You should now have access to the mongo shell via the `mongo` command and the daemon via the `mongod` command. Before we create the database and
   set up the daemon, we have to prepare a configuration file and decide where the database should be located. The location will be referred to
   as `<YourMongoPath>`, and you put the configuration file `mongod.conf` there. It should look somewhat similar to the yaml file below, but you can
   also tune it according to the options given [here](https://docs.mongodb.com/v4.0/reference/configuration-options/). Important settings
   are `storage.dbPath` (where your database will be located), `processManagement.fork` (If true, the process will run in the background and not block
   your shell), `net.port` (use the standard port for MongoDB: 27017), and `security.authorization` (should be disabled in the beginning to set up the
   users). Note that you must use spaces and not tabs in YAML!

Example mongod.conf file:<a name="mongodconf"></a>

 ~~~
storage:  
    dbPath: "<YourMongoPath>/data/db"  
    journal:  
        enabled: true  
        commitIntervalMs: 100  
    directoryPerDB: true  
    engine: "wiredTiger"  
    wiredTiger:  
        engineConfig:  
           cacheSizeGB: 4.0  
           journalCompressor: "snappy"  
           directoryForIndexes: false  
systemLog:
    verbosity: 0  
    quiet: false  
    traceAllExceptions: false  
    path: "<YourMongoPath>/mongod.log"  
    logAppend: true  
    logRotate: "reopen"  
    destination: "file"  
    component:  
        accessControl:  
            verbosity: 0  
        command:  
            verbosity: 0
 processManagement:  
    fork: true  
    pidFilePath: "<YourMongoPath>/pidfile"
net:  
    port: 27017  
    bindIp: localhost # Can also use a comma separated list  
    bindIpAll: false # If true can access from everywhere.  
    ipv6: true  
security:  
    authorization: disabled
~~~

3. Now, create the folder to hold the database files by typing:  
   `cd <YourMongoPath>; mkdir data; mkdir data/db` where you can choose `<YourMongoPath>` however you like, and then start the mongo daemon using the
   configuration file by executing: `mongod --config mongod.conf` which should result in something ending
   with : `child process started successfully, parent exiting` and produce a `mongod.log` and a `pidfile`. If you have some problems, check the
   logfile.

4. Now, we will use the mongo shell to create two users for the databases. We start by going in the shell and switching to the admin database (which
   we will use for authorization):  
   `mongo`
   `use admin`
   Now we will create an adminUser and a readOnlyUser by typing:
   `db.createUser({user: "<RootUser>", pwd: "<RootPassword>", roles: [{role: "root", db: "admin" }]})`
   and `db.createUser({user: "<ReadUser>", pwd: "<ReadPassword>", roles: [{role:"readAnyDatabase", db: "admin" }]})`
   Now exit the mongo shell:
   `exit`

5. Stop the mongo daemon by executing `mongod --shutdown --dbpath <YoutMongoPath>/data/db` and activate authorization in the `mongod.conf` file by
   changing the last line to `authorization: enabled`. Now, your database is password-protected, and you have to the username and password you created
   before in the mongo shell. To do that, while in the mongo shell, switch to the admin database by executing `use admin` and then authenticate
   with `db.auth(<RootUser>, <RootPassword>)`. Then exit the mongo shell again.

6. It might be a good idea to define some aliases in your `.bashrc` or `.myaliases` files to start and stop the mongod daemon, so you can do it
   quickly if needed. E.g.:

```
alias mongo_start="mongod -f <YourMongoPath>/mongod.conf"
alias mongo_stop="mongod --shutdown --dbpath <YourMongoPath>/data/db"
```

[Back to top](#toc)

### Installing surfflow<a name="surfflowinstall"></a>

1. You have two options here. For both options, first make sure that you are in the `surfflow` environment you created. If you want to install the
   latest published version of surfflow, you can simply run `pip install surfflow`. If you want to install the latest development version however, you can
   clone the repository by typing `git clone https://github.com/fyalcin/surfflow.git` and then run `pip install -e .` in the root directory of the
   repository. This will install the package in editable mode, so you can make changes to the code and see the changes immediately.
2. Select a location where you want to have your surfflow files located. We will assume this is `<YourPath>`. Now create two
   subfolders: `<YourPath>/config` for your configuration files and `<YourPath>/pps` for your pseudopotentials.

[Back to top](#toc)

### Configuration Files<a name="configfiles"></a>

surfflow comes with a nifty command line tool to make the generation of configuration files easy for the end user. After you have installed `surfflow`,
you can use the command `surfflow_tool config init` to start the configuration wizard that will run you through the process of defining the
configuration files and setting necessary `conda` environment variables. This tool creates a directory depending on user input (default is
at `~/config`), and all the necessary configuration files will be placed in this folder. This folder from now on will be referred to
as `<ConfigPath>`.

If you already have configuration files for `fireworks` in a directory, you can still use this tool to set the necessary conda environment variables
by pointing to this folder and skipping the config file generation part of the wizard. You can also use the command `surfflow_tool config init --help`
to see the options for the configuration wizard.

[Back to top](#toc)

### Configuring FireWorks for Job Schedulers<a name="configurefw"></a>

If the computer or cluster where you are installing surfflow has a job scheduler, you have to set up a file called `my_qadapter.yaml` in
the `<ConfigPath>` directory. `<SchedulerType>` can be PBS/Torque, SLURM, SGE, or IBM LoadLeveler. `pre_rocket` and `post_rocket` are optional
commands to be run in the job script before and after the workflow is executed; this is e.g. useful for loading and unloading modules. You probably
will have to activate your conda environment here. The `--timeout` option tells the job to stop pulling new FireWorks from the Launchpad after `<sec>`
number of seconds, which should of course be smaller than the walltime. E.g. if you have an allowed walltime of 72 hours, set the `--timeout` option
to e.g. 172800 (2 days in seconds).  
It is important to note that this type of rocket_launch will only work if the compute nodes on your cluster can access the MongoDB database, which is
a problem for many clusters since only the login nodes have access to the internet and firewall rules are strict. Possible solutions are
described [here](https://materialsproject.github.io/fireworks/offline_tutorial.html). The `my_qadapter.yaml` file might look something like this:

~~~
_fw_name: CommonAdapter  
_fw_q_type: <SchedulerType>
rocket_launch: rlaunch -c <YourPath>/config rapidfire --timeout <sec>
nodes: <NrOfNodes>
walltime: <hh:mm:ss>
queue: null  
account: null  
job_name: null  
pre_rocket: <Custom commands to load modules and conda etc.>
post_rocket: <Custom commands to unload modules etc.>
logdir: <YourPath>/logs
~~~

[Back to top](#toc)

### Configuring pymatgen<a name="configurepymatgen"></a>

When running VASP calculations FireWorks relies heavily on   [pymatgen](https://pymatgen.org/)
and [Custodian](https://materialsproject.github.io/custodian/). Some configuration of pymatgen is required:

1. We assume that you have a folder `<EXTRACTED_VASP_POTCAR>` where you have the VASP pseudopotentials (delivered with the VASP package) already
   extracted. Type `pmg config -p <EXTRACTED_VASP_POTCAR> <YourPath>/pps` to put these potentials in an order where pymatgen can find them. The final
   file structure should look something like this (you maybe have to rename the directories in the pps folder):

~~~
pps
├── POT_GGA_PAW_PBE  
│   ├── POTCAR.Ac.gz  
│   ├── POTCAR.Ac_s.gz  
│   ├── POTCAR.Ag.gz  
│   └── ...  
├── POT_GGA_PAW_PW91  
│   ├── POTCAR.Ac.gz  
│   ├── POTCAR.Ac_s.gz  
│   ├── POTCAR.Ag.gz  
│   └── ...  
└── POT_LDA_PAW  
    ├── POTCAR.Ac.gz  
    ├── POTCAR.Ac_s.gz  
    ├── POTCAR.Ag.gz  
    └── ...
~~~

2. Now, we have to set a config variable (it will be a file .`pmgrc.yaml` in your home folder) so pymatgen can find the potentials and add your
   default functional as well (this could also be PBE_54) if you have this potential family and did not rename the folders in the previous step):  
   `pmg config --add PMG_VASP_PSP_DIR <YourPath>/pps`
   `pmg config --add PMG_DEFAULT_FUNCTIONAL PBE`
3. For the integration of the Materials Project REST API, you should register for free at the
   website [https://materialsproject.org/](https://materialsproject.org/) and get an API Key from your dashboard there. Put it in the configuration
   file:  
   `pmg config --add PMG_MAPI_KEY <Your_API_Key>`

[Back to top](#toc)

## Testing your installation of FireWorks<a name="testing"></a>

- Try to load a workflow (for a simple structure relaxation of Si) to the launchpad with `atwf add -l vasp -p wf_structure_optimization -m mp-149` and
  check that it has been added with `lpad get_wflows`. This should result in something like:

~~~
[  
 {  
  "state": "READY",  
  "name": "Si—1",  
  "created_on": "2020-02-27T14:44:42.634000",  
  "states_list": "REA"  
 },  
]
~~~

- Navigate to a scratch directory and run the workflow (without a scheduler) with `rlaunch rapidfire`
- Afterwards, you can run `lpad get_wflows` again to see that the state has changed from "READY" to "COMPLETED"
- It would probably be a good idea to continue with the [tutorial of Atomate](https://atomate.org/running_workflows.html) if you are not familiar with
  it already, but you can also jump straight into surfflow in the next section.

[Back to top](#toc)

## A note about atomate<a name="customatomate"></a>

Atomate supplies a bunch of Fireworks and workflows that are used in surfflow. However, we have added a little more custom functionality to it. Mainly
enabling the copy of the `vdw_kernel.bindat` file for `OptimizeFW`, `StaticFW` and `TransmuterFW` completely analogous to the way it was implemented
in `ScanOptimizeFW` already. Another small modification concerns the copying of `WAVECAR` files from a relaxation to a `StaticFW`, which was not
possible. surfflow automatically installs atomate from a forked repo, so those changes are implemented automatically during installation. Of course,
this is not an ideal solution, and we strive to have something similar implemented in the official atomate version, but so far, we were not
successful.

[Back to top](#toc)

## Using surfflow<a name="using"></a>

### Running a workflow<a name="runningawf"></a>

Main workflow is located in the `surfflow.workflows.main` module. One can submit workflows by either defining some input parameters in a script and
calling the workflow from this module, or by using the same command line tool we used for the configuration files. The latter is the recommended way
but we still provide a short example of the former in the `surfflow.tests.surface_energy_test` module. The command for submitting workflows
is `surfflow_tool submit` and has the following options:

```
  -h, --help            show this help message and exit
  -m MPID, --mpid MPID  MPID of the material to submit.
  --miller MILLER [MILLER ...]
                        Miller indices of the surface to submit.
  -f FILENAME, --file FILENAME
                        Path to the file containing the input parameters to use.
  -sg SG_PARAMS, --sg-params SG_PARAMS
                        SlabGenerator parameters given in a dict format.
```

The `--mpid` option is the only mandatory one. It is the Materials Project ID of the material to submit. The `--miller` option is a list of Miller
indices of the surface to submit. If not given, all surfaces up to a `max_index` of the material will be submitted. The `--file` option is the path to
a file containing the input parameters to use. The `--sg-params` option is a dictionary of parameters that define the slab geometries.

We provide an example of an input file that can be used with the `--file` option in `surfflow/tests/inp.yaml`, which you can submit with the
command `surfflow_tool submit --file <FullPathToInputFile>`.

For all the options, the inputs are checked against default values defined in `surfflow/defaults.json` and set accordingly if not supplied by the user.
If one wants to override the defaults given, the recommended way is to create a file `user_defaults.json` in `<ConfigPath>` and change the values
there. User may refer to the `surfflow/defaults.json` for the correct format for defining these parameters.

### Explanation of input parameters<a name="inputparams"></a>

Here, we will go over the input parameters in the `surfflow/defaults.json` file

- `comp_params` - Parameters used in VASP calculations
    - `use_vdw`: Whether to use the VdW-DF functional. If `true`, the `vdw_kernel.bindat` file will be copied to the working directory of the
      calculation. Default: `false`
    - `functional`: The functional to use. Default: `PBE`
    - `use_spin`: Whether to use spin-polarized calculations. Default: `true`
    - `encut`: The cutoff energy in eV. Default: `520`
    - `k_dens`: The k-point density in Å<sup>-1</sup>. Default: `6`
    - `is_metal`: Whether the material is a metal. Default: `null`
    - `epsilon`: The dielectric constant of the material. Default: `null`


- `sg_params` - Parameters used in surface generation
    - `symmetrize`: Whether to symmetrize the slab. Default: `false`
    - `slab_thick`: The thickness of the slab in number of layers. Default: `10`
    - `vac_thick`: The thickness of the vacuum in Å. Default: `20`
    - `primitive`: Whether to perform a primitive cell transformationon the slabs. Default: `true`
    - `lll_reduce`: Whether to perform a LLL reduction on the slabs. Default: `false`
    - `tol`: The tolerance for determining the layers in Å. Default: `0.1`
    - `max_index`: The maximum Miller index to consider if `miller` isn't supplied. Default: `3`
    - `max_normal_search`: The maximum number of normal vectors to search for the slab. Default: `max` (highest absolute value of individual Miller
      indices)
    - `resize`: Whether to resize the slabs to the desired thickness. Default: `true`
    - `preserve_terminations`: Whether to preserve the terminations of the slabs while resizing. Default: `true`
    - `min_thick_A`: The minimum thickness of slabs in Å. Default: `10`
    - `minimize_structures`: Whether to try all combinations of `primitive`, `lll_reduce`, and `max_normal_search` to find the smallest slabs.
      Default: `true`
    - `match_ouc_lattice`: Whether to match the lattice of the slab to the lattice of the oriented bulk structure. Default: `true`
    - `calculate_bonds`: Whether to calculate the bonds broken while generating the slabs. Default: `true`
    - `max_nsites`: The maximum number of sites in the slabs. Default: `100`
    - `nn_method`: The method to use for finding the nearest neighbors. Default: `all`
    - `weight`: The weight of the bonds broken while generating the slabs. Default: `BVS`


- `sg_filter`: Dictionary of options to filter the slabs.
    - `method`: Method to use for filtering the slabs based on the bonds broken. Default: `all`, which means nothing is filtered out.
      Available methods are `bvs_threshold`, `bvs_min_N`, and `bvs_min_N_hkl`.
        - `bvs_threshold`: Percentage based filtering. Only slabs with BVS values within a given percentage of the minimum BVS slab are kept.
        - `bvs_min_N`: Number based filtering. Slabs are ordered by BVS values and the slabs with the lowest `N` BVS values are kept.
        - `bvs_min_N_hkl`: Number based filtering within each Miller index. Same thing as `bvs_min_N` but the slabs are grouped by Miller index first
          and the slabs with the lowest `N` BVS values are kept for each group.
    - `bvs_param`: Parameter `N` to use for all the methods. Should be given as a fraction for `bvs_threshold` and an integer for the other methods.


- `bulk_coll_suffix`: The suffix to use for the bulk collection. Default: `bulk_data`
- `surfflow_coll_suffix`: The suffix to use for the surfflow collection. Default: `surfflow_data`
- `db_file`: The path to the database file. Default: `auto` (uses the `db_file` in the `FW_CONFIG_FILE`)
- `high_level`: Name of the high_level database to save the results to. Default: `true` (uses the `high_level` in the `FW_CONFIG_FILE`)
- `functional`: The functional to use for the high-level calculations. Default: `PBE` (this is repeated but it is required for now for the workflow to
  work)
- `fake_vasp`: Whether to use fake VASP. Instead of running VASP, it copies inputs and outputs from a directory in the module path. Will give nonsense
  results for surface energies. Commonly used for testing purposes. Use this if you haven't set up VASP yet but would like to test the workflow.
  Default: `false`

[Back to top](#toc)

### Looking at the results<a name="results"></a>

The results of surfflow are saved in a separate MongoDB database, which nevertheless is hosted on the same server (note that the `directoryPerDB: true`
line in [`mongod.conf`](#mongodconf) assures that the results are stored in a different folder). This database of results is set in
the [db.json](#dbjson) file under the key "high_level", while the "database" database is used by FireWorks to store all kinds of stuff. The results
are stored in different [collections](https://docs.mongodb.com/manual/core/databases-and-collections/#databases), separated for the functional used (
PBE, or SCAN) and then split between bulk, slab, and interfaces results. Data in the surfflow database can be queried, of course, directly from
the [mongo shell](https://docs.mongodb.com/manual/mongo/), but it is probably more useful to use the python interface to
MongoDB, [pymongo](https://pymongo.readthedocs.io/en/stable/). Some functions to aid with this are provided in the `surfflow.utils` module. To look
quickly at single results or get a feel for how the data is structured in the database, it might be beneficial to install a GUI for MongoDB,
like [Compass](https://www.mongodb.com/products/compass). (Note that you can use the web-GUI of Atlas when you are using this cloud-based solution
instead of a local installation of MongoDB.)

Also, note that there is a web-GUI provided by FireWorks where you can check out the Workflows and Fireworks in your FireWorks database. Just
type `lpad web_gui` for that.

[Back to top](#toc)

