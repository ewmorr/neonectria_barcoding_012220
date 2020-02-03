install for 3prime trim mod of itsxpress

First install a qiime2 conda env (instructions here: https://docs.qiime2.org/2019.10/install/native/)
```
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-osx-conda.yml
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-osx-conda.yml
rm qiime2-2019.10-py36-osx-conda.yml
conda activate itsxpress_3p
```
install dependencies
```
conda install hmmer bbmap vsearch
```
clone itsxpress from ewmorr github,  switch to 3p_trim branch, and then install
```
git clone https://github.com/ewmorr/itsxpress.git
cd ~repo/itsxpress
git branch
git checkout 3p_trim
pip install -e .
```
