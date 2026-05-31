We are using a custom container that includes all R and Python packages required for the visualization modules. The container is built from the `groupD.def` definition file.

The resulting `.sif` image is too large to be stored in the GitHub repository, so it is not included.

To obtain the container, you can either:

* Copy `groupD.sif` from the `containers/` directory (if available).
* Retrieve it from the cluster.
* Build it locally using:

```bash
sudo apptainer build groupD.sif groupD.def
```
**Note:** Building the container may take several minutes, as a number of R and Python dependencies need to be installed during the process.
