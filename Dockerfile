FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    libxml2-dev \
    libz-dev \
    libglib2.0-dev \
    libffi-dev \
    python3-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . /app

EXPOSE 8000

CMD ["python", "server.py"]

