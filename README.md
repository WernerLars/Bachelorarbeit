###Installation for Windows 10

#First Install following Software
Matlab 2020b (needed for IronClust,Waveclus)
	- Signal Toolbox
	- Image Processing Toolbox
	
Visual Studio 19: https://visualstudio.microsoft.com/de/downloads/
	- Workloads choose Desktop C++
		- MSVC v142 VS 2019 C++ x64 Buildtools
		- Windows 10 SDK (10.0.183.62.0)
		- Just in Time Debugger
		- C++ Profile Tools
		- C++ CMake Tools for WINDOWS
		- C++ ATL v142 Buildtools

MPICH:  https://www.microsoft.com/en-us/download/details.aspx?id=100593	 	(needed for Spyking Circus)
	you need do install msmpisetup.exe and msmpisdk.msi

# You will need to install it on the System Drive (C:\) (Spyking Circus has Path ValueErrors if you don't)

#Open CMD in the cloned Gitlab Repository

#Create venv
python -m venv venv

#Start venv
venv\Scripts\activate.bat

#Run following code
python -m pip install --upgrade pip
pip install -r Requirements.txt

#HDSCAN ON WINDOWS HAS A BUG WITH NUMPY (MIGHT GET FIXED IN THE FUTURE!!) SO DO THE FOLLOWING

git clone https://github.com/scikit-learn-contrib/hdbscan
cd hdbscan
git checkout ea81c103a49f47fa9388544038e85fb833ff1233

#OPEN IN WINDOWS EXPLORER hdbscan\pyproject.toml with an Editor and edit numpy version to 1.19.3:

[build-system]
requires = [
  "setuptools",
  "wheel",
  "cython",
  "numpy==1.19.3"
]

#In Console
python setup.py install
cd..

#Now install the Sorter:

pip install -r Software.txt
git clone https://github.com/flatironinstitute/spikeforest_recordings
cd spikeforest_recordings
git checkout dc636f392c0a870ac8b69ac7a7ed978161df39fd					# For the Commit Version I have tested the code
cd..
git clone https://github.com/flatironinstitute/ironclust
cd ironclust
git checkout e46c00e849a822a6d432fca334d9b89ad7b5da4b					# For the Commit Version I have tested the code
cd..
git clone https://github.com/csn-le/wave_clus
cd wave_clus
git checkout bfa34eaabeab75429d6e233b58bd1583c1872205					# For the Commit Version I have tested the code

#Now Close CMD

#Tridesclous installs older version of pyqtgraph which has an old version of ptime that is incompatilbe with Python Version 3.8.6 ( time has no attribute clock )
#Fix: Replace ptime.py in the Folder venv\Lib\site-packages\pyqtgraph with fixed version in my repository

#Set Systemvariables for:

KACHERY_STORAGE_DIR 	# InstalledPath/data	
IRONCLUST_PATH		# InstalledPath/ironclust
WAVECLUS_PATH		# InstalledPath/wave_clus

###Start Program

#Open CMD

#Start venv
venv\Scripts\activate.bat

#Run Jupyter Notebook
jupyter notebook

#Choose Spike_Sorting_Pipeline.ipynb