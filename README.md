# photozpy
The photo-z analysis pipeline with python

## Installation

1. It's recommended to use photozpy in a conda environment. Follow the instructions on the [Anaconda website](https://docs.anaconda.com/free/anaconda/install/index.html) for installation.

2. Create conda environment
	Open your terminal
	```bash
	conda create --name photozpy pip python=3.9.19
	```

3. Activate conda environment
	In the terminal, run the following command to activate the photozpy environment
	```bash
	conda activate photozpy
	```

4. Clone the photozpy GitHub to your local machine
	In the terminal, run the following command to download the photozpy
	```bash
	git clone https://github.com/Yong2Sheng/photozpy.git
	```
	Note that you might need to setup your ssh keys if it's your first time to clone a remote repository to your local machine.

5. Enter the *photozy* folder and install `photozpy` module into the photozpy environment
	```shell
	>>> cd photozpy
	>>> pip install -e .
	```

6. Once the installation is done, you can import and use `photozpy` module. Please go to the *docs* folder to try out the tutorial. 
