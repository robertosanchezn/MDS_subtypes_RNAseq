#!/bin/sh
#SBATCH --job-name=rstudio-server
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --output=rstudio.log

# Script based on https://rocker-project.org/use/singularity.html
# Changes to: 
# - use apptainer instead of singularity
# - run r and libraries inside the container

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(mktemp -d)

# Set R_LIBS_USER to an existing path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<"END"
#!/bin/sh
export R_LIBS_USER=/opt/conda/lib/R/library
# mkdir -p "${R_LIBS_USER}"
## custom Rprofile & Renviron (default is $HOME/.Rprofile and $HOME/.Renviron)
# export R_PROFILE_USER=/path/to/Rprofile
# export R_ENVIRON_USER=/path/to/Renviron
exec /usr/lib/rstudio-server/bin/rsession "${@}"
END

chmod +x ${workdir}/rsession.sh

r_path="/opt/conda/bin/R"
mkdir -p ${workdir}/rstudio

cat > ${workdir}/rstudio/rserver.conf <<END
rsession-which-r=${r_path}
END

export APPTAINER_BIND="${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/rstudio/rserver.conf:/etc/rstudio/rserver.conf"
# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0

export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD="m"
# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & apptainer commands
readonly PORT=$(python3 -c 'import socket; s=socket.socket(); s.bind(("", 58733)); print(s.getsockname()[1]); s.close()' || python3 -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:${HOSTNAME}:${PORT} ${APPTAINERENV_USER}@hpclogin.unav.es

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END

apptainer exec \
   --cleanenv \
   --scratch /run,/tmp,/var/lib/rstudio-server \
   --workdir ${workdir} \
   Apptainer.sif \
      rserver \
      --www-port ${PORT} \
      --auth-none=1 \
      --auth-pam-helper-path=pam-helper \
      --auth-stay-signed-in-days=30 \
      --auth-timeout-minutes=0 \
      --server-user=$(whoami) \
      --rsession-path=/etc/rstudio/rsession.sh

printf 'rserver exited' 1>&2
