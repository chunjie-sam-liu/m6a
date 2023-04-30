# Test CIMS on iclip dataset

https://zhanglab.c2b2.columbia.edu/index.php/ICLIP_data_analysis_using_CTK

## Install CTK

```bash
myenv='ctk'
conda create --yes --name $myenv

conda activate $myenv

conda config --env --append channels conda-forge
conda config --env --append channels bioconda

conda install --yes -c chaolinzhanglab ctk
```