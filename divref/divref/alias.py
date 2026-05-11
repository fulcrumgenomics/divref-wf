from typing import TypeAlias

HailPath: TypeAlias = str
"""Type alias for filesystem paths accepted by Hail: local, S3 (s3a://), GCS (gs://), or HDFS (hdfs://)."""  # noqa: E501
