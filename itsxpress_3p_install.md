install for 3prime trim mod of itsxpress

### stand-alone install (use when root priviledges are unavailable)

```
conda create -n trim_3p pip hmmer bbmap vsearch biopython
conda activate trim_3p
```
clone itsxpress from ewmorr github,  switch to 3p_trim branch, and then install
```
git clone https://github.com/ewmorr/itsxpress.git
cd ~/repo/itsxpress
git branch
git branch 3p_trim origin/3p_trim
git checkout 3p_trim
pip install --user -e .
```
`pip install --user` puts executable in `~/.local/bin` which may not be in $PATH by default
```
export PATH=~/.local/bin:$PATH
```

### To install within qiime2 env (requires root priviledges)
First install a qiime2 conda env (instructions here: https://docs.qiime2.org/2019.10/install/native/)
```
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-osx-conda.yml
conda env create -n itsxpress_3p --file qiime2-2019.10-py36-osx-conda.yml
```
linux
```
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml
conda env create -n itsxpress_3p --file qiime2-2019.10-py36-linux-conda.yml
```
cleanup and activate
```
rm qiime2-2019.10-py36-osx-conda.yml
rm qiime2-2019.10-py36-linux-conda.yml
conda activate itsxpress_3p
```
install dependencies (pip is needed on linux server so will not use system pip install)
```
conda install pip hmmer bbmap vsearch biopython
```
clone itsxpress from ewmorr github,  switch to 3p_trim branch, and then install
```
git clone https://github.com/ewmorr/itsxpress.git
cd ~/repo/itsxpress
git branch
git branch 3p_trim origin/3p_trim
git checkout 3p_trim
pip install --user -e .
```

