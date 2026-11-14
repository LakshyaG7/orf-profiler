FROM python:3.12-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir .

# Mount your data at /data, for example:
#   docker run --rm -v "$PWD/data:/data" orf-profiler scan /data/transcripts.fa -o /data/results
ENTRYPOINT ["orf-profiler"]
CMD ["--help"]
