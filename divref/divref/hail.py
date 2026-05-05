import os
from pathlib import Path

import hail as hl
import pyspark


def hail_init(
    gcs_credentials_path: Path, spark_driver_memory_gb: int = 1, spark_executor_memory_gb: int = 1
) -> None:
    """
    Initialize Hail with the S3A and GCS connectors and credential configuration.

    Sets ``GOOGLE_APPLICATION_CREDENTIALS`` so the JVM subprocess inherits it,
    then starts Hail with both the S3A connector JARs (``hadoop-aws`` and
    ``aws-java-sdk-bundle``) and the GCS connector JAR on the Spark classpath.
    AWS credentials are resolved via the standard
    ``DefaultAWSCredentialsProviderChain`` (env vars, ``~/.aws/credentials``,
    or IAM role).

    Args:
        gcs_credentials_path: Absolute path to a GCP Application Default Credentials
            JSON file. If the file exists and ``GOOGLE_APPLICATION_CREDENTIALS`` is not
            already set, it is exported to the environment before Hail starts. Required
            only if ``gs://`` paths are read; ignored otherwise.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark executor memory must be at least 1GB. Saw {spark_executor_memory_gb}GB."
        )

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )

    if gcs_credentials_path.exists() and "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ:
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(gcs_credentials_path)

    jars_dir = os.path.join(pyspark.__path__[0], "jars")
    gcs_jar = os.path.join(jars_dir, "gcs-connector.jar")
    hadoop_aws_jar = os.path.join(jars_dir, "hadoop-aws.jar")
    aws_sdk_bundle_jar = os.path.join(jars_dir, "aws-java-sdk-bundle.jar")

    if not os.path.exists(gcs_jar):
        raise FileNotFoundError(
            f"GCS connector JAR not found at {gcs_jar}. Run 'pixi run setup-gcs' to download it."
        )
    for jar in (hadoop_aws_jar, aws_sdk_bundle_jar):
        if not os.path.exists(jar):
            raise FileNotFoundError(
                f"S3 connector JAR not found at {jar}. Run 'pixi run setup-s3' to download it."
            )

    cloud_jars = ",".join([gcs_jar, hadoop_aws_jar, aws_sdk_bundle_jar])

    hl.init(
        spark_conf={
            "spark.jars": cloud_jars,
            "spark.driver.extraClassPath": cloud_jars,
            "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
            "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",  # noqa: E501
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",  # noqa: E501
            # Random-read optimizations for Hail/Parquet workloads. S3A defaults to
            # sequential fadvise, which re-opens the HTTP connection on every backward
            # seek; Hail does many seeks per partition, so `random` is much faster.
            "spark.hadoop.fs.s3a.experimental.input.fadvise": "random",
            "spark.hadoop.fs.s3a.readahead.range": "64K",
            # Defaults (15 / 10) bottleneck partition-parallel reads.
            "spark.hadoop.fs.s3a.connection.maximum": "200",
            "spark.hadoop.fs.s3a.threads.max": "64",
        }
    )
