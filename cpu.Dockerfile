# cpu.Dockerfile

# 1) Start from a slim Python base
FROM python:3.10-slim

WORKDIR /opt/program

# 2) Copy & prune out any GPU‐specific pins from your requirements
COPY requirements-gpu.txt ./
RUN grep -Ev '^(torch|torchvision|torchaudio)' requirements-gpu.txt \
    > requirements-cpu.txt

# 3) Install system tools needed to build some Python packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc libpq-dev \
 && rm -rf /var/lib/apt/lists/*

# 4) Install CPU‐only PyTorch, TorchVision, and Torchaudio
RUN pip install --no-cache-dir \
      torch torchvision torchaudio \
      -f https://download.pytorch.org/whl/cpu/torch_stable.html

# 5) Install your other Python deps + Gunicorn
RUN pip install --no-cache-dir \
      -r requirements-cpu.txt \
      gunicorn

# 6) Copy in your code
COPY . .

# 7) Create the same “serve” entrypoint as your GPU image
RUN printf '#!/bin/bash\n\
if [ "$1" = "serve" ]; then\n\
  echo "[$(date)] Starting inference server on port 8080…"\n\
  exec gunicorn --bind 0.0.0.0:8080 run:app\n\
else\n\
  exec "$@"\n\
fi\n' \
  > /usr/local/bin/serve \
  && chmod +x /usr/local/bin/serve

ENTRYPOINT ["serve"]
CMD []
