# Running the OpenQuake Engine on multiple nodes (cluster/zmq)

This configuration is supported only by Linux.

## Overall architecture
The nodes must all be able to communicate with the OpenQuake Engine *DbServer*.
Both services run on a single "master" node.
Running OpenQuake on an *MPI cluster* is currently not supported. See the [FAQ](../faq.md#mpi-support) for more information.

## Initial install

### Pre-requisites

Have read [Installing on RedHat and derivatives](rhel.md) or [Installing on Ubuntu](ubuntu.md) (depending on the operating system been used).

### Master node
The `python3-oq-engine-master` package must be installed on the **master** node.

On **RHEL/CentOS** [EPEL](https://fedoraproject.org/wiki/EPEL) repository *must be configured and enabled* in the system.

### Worker nodes
On **worker** nodes  `python3-oq-engine-worker` must be installed **instead**.

## OpenQuake Engine 'master' node configuration File

### Enable zmq distribution

The following file (on all nodes) should be modified to enable
*zmq* support:

`/etc/openquake/openquake.cfg`

```
[distribution]
# enable celery only if you have a cluster
oq_distribute = zmq

[dbserver]
multi_user = true
file = /var/lib/openquake/db.sqlite3
# daemon bind address; must be a valid IP address
listen = 0.0.0.0
# address of the dbserver; can be an hostname too
# on multi-node cluster it must be the IP or hostname
# of the master node (on the master node cfg too)
host = w.x.y.z
port = 1907
receiver_ports = 1912-1920
authkey = somethingstronger

[zworkers]
host_cores = 192.168.2.1 -1, 192.168.2.2 -1, 192.168.2.3 -1, 192.168.2.4 -1, 192.168.2.5 -1
ctrl_port = 1909
```

Notice that the -1 in `192.168.2.1 -1` means that all the cores in
the worker with IP `192.168.2.1` will be used. You can use a number
between 0 and the maximum number of available core to limit the
resource usage. The engine will automatically start and stop zmq
processes on the worker nodes at each new calculation, *provided the
user openquake has ssh access to the workers*.

NB: when using the zmq mechanism you should not touch the parameter
`serialize_jobs` and keep it at its default value of `true`.


### Configuring daemons

The required daemons are:

#### Master node
- OpenQuake Engine DbServer
- OpenQuake Engine WebUI (optional)

### Monitoring zmq

`oq workers status` can be used to check the status of the worker nodes and the task distribution. An output like this is produced:

```
$ oq workers status
[('192.168.2.1', 1, 64), ('192.168.2.2', 7, 64), ('192.168.2.3', 7, 64)]
```

For each worker in the cluster you can see its IP and the cores which are
currently running with respect to the number of cores available (for instance
on the host 192.168.2.1 only 1 core of 64 is running, while in the other
two workers 7 cores are running each).

## Shared filesystem

OpenQuake 2.4 introduced the concept of _shared directory_ (aka _shared_dir_). This _shared directory_ allows the workers to read directly from the master's filesystem, thus increasing scalability and performance; starting with OpenQuake 3.3 this feature is **mandatory** on a multi-node deployment.

The _shared directory_ must be exported from the master node to the workers via a _POSIX_ compliant filesystem (**NFSv4** is recommended). The export may be (and _should_ be) exported and/or mounted as **read-only** by the workers.

As soon as the shared export is in place, the `shared_dir` config parameter in `openquake.cfg` must be configured to point to the path of the exported dir on the _master_ node and to where the export is mounted on each _worker_ node.

```
[directory]
# the base directory containing the <user>/oqdata directories;
# if set, it should be on a shared filesystem; for instance in our
# cluster it is /home/openquake; if not set, the oqdata directories
# go into $HOME/oqdata, unless the user sets his own OQ_DATADIR variable
shared_dir = /home/openquake
```

When `shared_dir` is set, the `oqdata` folders will be stored under `$shared_dir/<user>/oqdata` instead of `/home/<user>/oqdata`. See the comment in the `openquake.cfg` for further information.
You need then to give `RWX` permission to the `shared_dir` on _master_ to the `openquake` group (which is usually created by packages) and add all the cluster users to the `openquake` group. For example:

```bash
$ mkdir /home/openquake
$ chown openquake.openquake /home/openquake
# Enable setgid on the tree
$ chmod 2770 /home/openquake
```

On the workers the _shared_dir_ should be mounted as the `openquake` user too, or access must be given to the user running `celeryd` (which is `openquake` by default in the official packages).

Another possibility would be exporting the entire `/home` to the workers: in such case `oqdata` will have the default path `/home/<user>/oqdata` and setgid is not required. Please note that the `openquake` user on workers still needs to get access to the `oqdata` content, so make sure that permission are properly set (`traverse` on the user home and `read` access to oqdata).

```
[directory]
shared_dir = /home
```

## Network and security considerations

The worker nodes should be isolated from the external network using either a dedicated internal network or a firewall.
Additionally, access to the DbServer ports should be limited (again by internal LAN or firewall) so that external traffic is excluded.

The following ports must be open on the **master node**:

* 1907 for DbServer (or any other port allocated for the DbServer in the `openquake.cfg`)
* 1911 for the ZeroMQ streamer
* 1912-1920 for ZeroMQ receivers
* 8800 for the API/WebUI (optional)

The following port must be open on the **workers node**:

* 1909 for the ZeroMQ workerpools

The **worker nodes** must be able to connect to the master on port 1907.
Moreover the master must be able to access the workers via ssh.
This means that you have to generate and copy the ssh keys properly, and
the first time you must connect to the workers manually. Then the engine
will be able to start and stop zworker processes at each new calculation.

## Storage requirements

Storage requirements depend a lot on the type of calculations you want to run. On a worker node you will need just the space for the operating system, the logs and the OpenQuake installation: less than 20GB are usually enough. Workers can be also diskless (using iSCSI or NFS for example). Starting from OpenQuake 2.3 the software and all its libraries are located in `/opt/openquake`.

On the master node you will also need space for:
- the users' **home** directory (usually located under `/home`): it contains the calculations datastore (`hdf5` files located in the `oqdata` folder)
- the OpenQuake database (located under `/var/lib/openquake`): it contains only logs and metadata, the expected size is tens of megabyte
- the temporary folder (`/tmp`). A different temporary folder can be customized via the `openquake.cfg`

On large installations we strongly suggest to create separate partition for `/home` and `/tmp` (or any custom temporary folder set in the `openquake.cfg`.


## Swap partitions

Having swap active on resources dedicated to the OpenQuake Engine is _strongly discouraged_ because of the performance penality when it's being used and because how python allocates memory. In most cases (when memory throughput is relevant) is totally useless and it will just increase by many orders of magnitude the time required to complete the job (making the job actually stuck).


## Running calculations

Jobs can be submitted through the master node using the `oq engine` command line interface, the API or the WebUI if active. See the documentation about [how to run a calculation](../running/unix.md) or about how to use the [WebUI](../running/server.md)
