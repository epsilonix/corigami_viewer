#!/bin/sh
aws s3 sync --only-show-errors s3://corigami-data/ /app/corigami_data/
exec "$@"        # hand off to CMD or gunicorn