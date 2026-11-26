# Steps

When all tests done and successful and PR is approved by codeowners, follow these steps:

1. Ensure that you update the version in `microSALT/__init__.py`
1. Select "Squash and merge" to merge branch into default branch (master).
1. Deploy master to stage:
    1. Log in to appropriate server `ssh <server.scilifelab.se>`
    1. The deploy script must run as `hiseq.clinical`: `sudo -iu hiseq.clinical`, then `us`
    1. Install microSALT:
        ```shell
        bash /home/proj/stage/servers/resources/hasta.scilifelab.se/update-microsalt.sh -e stage -b master
        ```
        Flags: `-e` deployment environment, `-b` branch, `-a` run the Alembic migration (only when
        the release includes one), `-d` reinstall all dependencies.
1. Deploy master to production:
    1. Log in to appropriate server `ssh <server.scilifelab.se>`
    1. `sudo -iu hiseq.clinical`, then `up`
    1. Install microSALT:
        ```shell
        bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt.sh -e production -b master
        ```
