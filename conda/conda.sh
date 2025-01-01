conda create -n MGGNN -c conda-forge -c pytorch \
python=3.8 \
numpy=1.22.3 \
pandas=1.4.2 \
scipy=1.8.1 \
scikit-learn=1.1.1 \
scikit-misc=0.1.4 \
tqdm=4.64.0 \
pytorch">=1.8.0" \
cudatoolkit">=10.2" \
matplotlib=3.4.2 \
scanpy=1.9.1 \
anndata=0.8.0 \
r-mclust=5.4.10 \
rpy2 \
pot

conda activate MGGNN
python -m pip install spatialentropy
python -m pip install GraphST