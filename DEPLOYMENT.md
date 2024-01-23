# Steps

When all tests done and successful and PR is approved by codeowners, follow these steps:

1. Ensure that you update the version in `microSALT/__init__.py`
2. Select "Squash and merge" to merge branch into default branch (master).
3. Deploy master to stage:
   1. Log in to appropriate server `ssh <server.scilifelab.se>`
   2. `us`
   3. If you have updated the `environment.yml` or `requirements.txt`, run the following:
      ```shell
      conda env create -y -n S_microSALT -f https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/environment.yml
      conda activate S_microSALT
      pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/requirements.txt
      conda deactivate
      ```
   4. ```shell
      bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt-stage.sh master
      ```
   5. Make sure that the installation was successful
6. Deploy master to production
   1. Log in to appropriate server `ssh <server.scilifelab.se>`
   2. `up`
   3. If you have updated the `environment.yml` or `requirements.txt`, run the following:
      ```shell
      conda env create -y -f https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/environment.yml
      conda activate P_microSALT
      pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/requirements.txt
      conda deactivate
      ```
   4. ```shell
      bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt-prod.sh
      ```
   5. Make sure that the installation was successful