# Steps

When all tests done and successful and PR is approved by codeowners, follow these steps:

1. Ensure that you update the version in `microSALT/__init__.py`
2. Select "Squash and merge" to merge branch into default branch (master).
3. Deploy master to stage:
    1. Log in to appropriate server `ssh <server.scilifelab.se>`
    2. `us`
    3. Install microSALT using **one** of the following commands
        1. If you have updated the `environment.yml`:
            ```shell
            sudo -iu hiseq.clinical
            us
            bash /home/proj/production/servers/resources/hasta.scilifelab.se/install-microSALT-stage.sh master
            ```
        2. Otherwise run:
            ```shell
            bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-tool-stage.sh -e S_microSALT -t microSALT -b master
            ```
    4. Make sure that the installation was successful
4. Deploy master to production:
    1. Log in to appropriate server `ssh <server.scilifelab.se>`
    2. `up`
    3. Install microSALT using **one** of the following commands
        1. If you have updated the `environment.yml`:
            ```shell
            sudo -iu hiseq.clinical
            up
            bash /home/proj/production/servers/resources/hasta.scilifelab.se/install-microSALT-prod.sh
            ```
        2. Otherwise run:
            ```shell
            bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-tool-prod.sh -e P_microSALT -t microSALT -b master
            ```
    4. Make sure that the installation was successful
