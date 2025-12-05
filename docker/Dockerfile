FROM md_software:latest

RUN apt-get update && apt-get install -y \
    openbabel \
    build-essential python3-dev gfortran git wget curl \
    && rm -rf /var/lib/apt/lists/*

# Force Python to use unbuffered output
ENV PYTHONUNBUFFERED=1

# Download Conda

ARG CONDA_DIR=/opt/conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh
ENV PATH="$CONDA_DIR/bin:$PATH"

# Accept Conda's Terms and Conditions

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    conda clean -afy

ARG CONDA_ENV_NAME=dynagent
ENV CONDA_DEFAULT_ENV=$CONDA_ENV_NAME
ENV PATH="$CONDA_DIR/envs/$CONDA_ENV_NAME/bin:$PATH"


# Create a new user for the agent
ARG AGENT_USER=beautifulagent
RUN useradd -m $AGENT_USER

WORKDIR /app
ENV SANDBOX=/app/sandbox
ENV AGENT_LOGS=/app/agent_logs

COPY . /app
RUN chmod -R a+rX,u-w,go-w /app

# Create a sandbox for the agent
RUN mkdir -p $SANDBOX && chown -R $AGENT_USER:$AGENT_USER $SANDBOX
RUN mkdir -p $AGENT_LOGS && chown -R $AGENT_USER:$AGENT_USER $AGENT_LOGS
RUN mkdir -p .venv && chown -R $AGENT_USER:$AGENT_USER .venv

COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["python", "/app/main.py"]