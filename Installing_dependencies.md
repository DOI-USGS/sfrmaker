## Install python dependencies
**Download and install the [Anaconda python distribution](https://www.anaconda.com/distribution/)**.
Open an Anaconda Command Prompt on Windows or a terminal window on OSX.
From the root folder for this repo (which contains `requirements.yml`), create a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) with the required python packages listed in `requirements.yml`:

```
conda env create -f requirements.yml
```

**if you want to contribute to sfrmaker**, use `ci/requirements.yml` instead, which includes additional packages for testing and documentation. Note that this environment is called `test` instead of `sfrmaker`.

### USGS users:
##### If you are on the USGS network, using Windows, and you get this error message:
```
CondaHTTPError: HTTP 500 INTERNAL ERROR for url <https://repo.anaconda.com/pkgs/msys2/win-64/m2w64-gettext-0.19.7-2.tar.bz2>
Elapsed: 00:30.647993

An HTTP error occurred when trying to retrieve this URL.
HTTP errors are often intermittent, and a simple retry will get you on your way.
```
Uncommenting the following line from `requirements.yml` should work:  
```
- msys2::m2w64-gettext
```

This tells conda to fetch `m2w64-gettext` from the `msys2` channel instead. Note that this is only a dependency on Windows,
so it needs to be commented out on other operating systems (normally it wouldn't need to be listed, but the above HTTP 500 error indicates that installation from the default source location failed.)

##### If you are on the USGS network and get an SSL error message 
(something similar to `SSL: CERTIFICATE_VERIFY_FAILED`), you need to configure `pip` package installer to use the USGS certificate.   

* On Windows, create the file `C:\Users\<your username>\AppData\Roaming\pip\pip.ini`.
* On OSX, create `/Users/<your username>/Library/Application Support/pip/pip.conf`

Include the following in this file:

    [global]
    cert = <path to DOI certificate file (e.g. DOIRootCA2.cer)>
    
Once this file is made, either  
	a) 	Update the `sfrmaker` conda environment with `conda env update -f requirements.yml` or  
	b) 	Remove and reinstall it:  
	
	conda env remove -n sfrmaker
	conda env create -f requirements.yml
	
Note that when you are off the USGS network, you may have to comment out the line in the above pip configuration file to get `pip` to work.

Once the `sfrmaker` conda enviornment is successfully made, activate it, so that the version of python in it will be at the top of the system path:

```
conda activate sfrmaker
```
##### If you are on the USGS network, using Windows, and encountering persistant issues with creating the conda environment
You may have better luck trying the install off the USGS network (e.g. at home). It is unclear whether these issues are from the USGS network, conda, or problems with the Windows versions of the package dependencies.


### Installing packages dependencies that are in active development (for easy updating)
The following python packages are being actively developed, and may need to be updated for sfrmaker to work properly. This can be accomplished by removing and rebuilding the `sfrmaker` conda environment, as above. Alternatively, the packages can be cloned from GitHub, and then linked to the `sfrmaker` conda environment so that updates to their source code are run automatically.

* [**flopy**](https://github.com/modflowpy/flopy)
* [**gis-utils**](https://github.com/aleaf/gis-utils)


#### Cloning a package
Urls to clone each package can be found at the links above, under the "Clone or Download" button.
Make sure you have [git installed](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). Then, in a command window pointed at the folder where you'd like to store the repository:

```
git clone https://github.com/modflowpy/flopy.git
```
(to clone flopy, for example)
#### Installing by linking the source code
Then, with the `sfrmaker` conda enviornment activated, `cd` into the repository folder and run the install script:

```
conda activate sfrmaker
cd flopy
pip install -e .
```
#### Incorporating updates
For example, to update flopy: from the flopy repository folder:

```
git fetch origin master
git rebase origin/master
```
This moves the state of your local copy to that of the remote (master branch) on code.usgs.gov, applying any changes you've made on top of the latest version on code.usgs.gov. If this leads to merge conflicts, you must resolve these before submitting your pull request. If you have uncommitted changes, you will need to git stash them prior to updating. This will effectively store your changes and they can be reapplied after updating.
